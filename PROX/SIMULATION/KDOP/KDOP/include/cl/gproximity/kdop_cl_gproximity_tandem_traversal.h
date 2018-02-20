#ifndef KDOP_CL_GPROXIMITY_TANDEM_TRAVERSAL_H
#define	KDOP_CL_GPROXIMITY_TANDEM_TRAVERSAL_H

#include <cl/kdop_cl_kernels_path.h>
#include <cl/kdop_cl_tandem_traversal.h>

#include <dikucl.hpp>
#include <dikucl_command_queue_manager.hpp>
#include <dikucl_context_manager.hpp>
#include <dikucl_kernel_manager.hpp>
#include <dikucl_util.h>

#include <contacts/geometry_contacts_callback.h>
#include <types/geometry_dop.h>

#include <kdop_tags.h>
#include <kdop_test_pair.h>
#include <kdop_tree.h>

#include <mesh_array.h>
#include <mesh_array_t4mesh.h>
#include <mesh_array_vertex_attribute.h>

#include <util_profiling.h>
#include <util_log.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdlib.h>
#include <string>
#include <utility>

#define GLOBAL_WORK_QUEUE_CAPACITY 50000

#ifdef __APPLE__
#define TT_THREADS 1
#else
#define TT_THREADS 64
#endif // __APPLE__
#define TT_WORK_GROUPS 128
#define TT_LOCAL_WORK_QUEUE_CAPACITY 600
#define TT_LOCAL_WORK_QUEUE_INIT_ITEMS TT_THREADS * 3
#define TT_GLOBAL_EXACT_TESTS_CAPACITY 16000000
#define TT_IDLE_THRESHOLD TT_WORK_GROUPS * 2 / 3

#ifdef __APPLE__
#define BW_THREADS 1
#else
#define BW_THREADS 512
#endif // __APPLE__
#define BW_WORK_GROUPS 1

#define ET_THREADS 64
#define ET_WORK_GROUPS 128

namespace kdop
{

  namespace details
  {

    namespace cl
    {

      template< typename KernelI >
      class gProximityKernelWorkItemGenerator
      : public KernelWorkItemGenerator<KernelI>
      {

      public:

        gProximityKernelWorkItemGenerator()
        {
        }

        void add_roots(std::vector<KernelI> a, std::vector<KernelI> b, KernelI tp, KernelI /* level */) {
          this->m_q.push(typename KernelWorkItemGenerator<KernelI>::RootPairs());
          this->m_q.back().a = std::vector<KernelI>(a);
          this->m_q.back().b = std::vector<KernelI>(b);
          this->m_q.back().tp = tp;
        }

        // 2016-06-22 Kenny code reivew:  'generate_kernel_work_items' hides
        // overloaded virtual function. Declared here: different number
        // of parameters (9 vs 7)
        size_t generate_kernel_work_items(
                                          size_t n,
                                          KernelWorkItem<KernelI> **out_kernel_work_items,
                                          size_t global_mem_cacheline_size,
                                          size_t /* max_bvtt_degree */,
                                          size_t /* max_bvtt_height */,
                                          size_t /* max_global_work_size */,
                                          size_t stride = 1
                                          )
        {
          if(this->m_q.empty()) {
            *out_kernel_work_items = NULL;
            return 0;
          }

          *out_kernel_work_items = new KernelWorkItem<KernelI>[n * stride];

          size_t actual = 0;
          while(actual < n && !(this->m_q.empty())) {
            for(; this->m_a < this->m_q.front().a.size() && actual < n; ++(this->m_a)) {
              for(; this->m_b < this->m_q.front().b.size() && actual < n; ++(this->m_b)) {
                (*out_kernel_work_items)[actual * stride].a = (KernelI) this->m_q.front().a[this->m_a];
                (*out_kernel_work_items)[actual * stride].b = (KernelI) this->m_q.front().b[this->m_b];
                (*out_kernel_work_items)[actual * stride].tp = (KernelI) this->m_q.front().tp;
                ++actual;
              }
              if(actual < n) {
                this->m_b = 0;
              } else {
                break;
              }
            }
            if(actual < n) {
              this->m_a = 0;
              this->m_q.pop();
            }
          }

          return actual;
        }



        // 2016-06-22 Kenny code reivew:  'cleanup_generated_work_items' hides
        // overloaded virtual function, declared here: different number of
        // parameters (2 vs 1)
        void cleanup_generated_work_items(
                                          KernelWorkItem<KernelI> * generated_kernel_work_items
                                          )
        {
          if(generated_kernel_work_items != NULL)
          {
            delete[] generated_kernel_work_items;
          }
        }

      };

    } // namespace cl

  } // namespace details

  template< typename V, size_t K, typename T, typename test_pair_container >
  inline void tandem_traversal(test_pair_container test_pairs,
                               dikucl::gproximity const & /* tag */,
                               size_t open_cl_platform = 0,
                               size_t open_cl_device = 0) {
    assert( ! test_pairs.empty() || !"tandem_traversal : test_pairs are empty" );

    cl_ulong tandem_traversal_kernel_time = 0;
    cl_ulong balance_work_kernel_time     = 0;
    cl_ulong exact_tests_kernel_time      = 0;

    size_t tandem_traversal_invocations = 0;
    size_t exact_test_invocations       = 0;

    cl_int err = CL_SUCCESS;

    ::dikucl::PlatformHandle *platform_handle = ::dikucl::PlatformManager::get_instance().get_platform(
                                                                                                       &err, open_cl_platform);
    CHECK_CL_ERR(err);

    ::dikucl::DeviceHandle *device_handle = ::dikucl::DeviceManager::get_instance().get_device(
                                                                                               &err, platform_handle, open_cl_device);
    CHECK_CL_ERR(err);
    size_t global_mem_cacheline_size =
    (size_t) device_handle->device.getInfo<CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE>(&err);
    CHECK_CL_ERR(err);

    ::dikucl::ContextHandle *context_handle = ::dikucl::ContextManager::get_instance().get_context(&err, device_handle);
    CHECK_CL_ERR(err);
    ::cl::Context context = context_handle->context;

    ::dikucl::CommandQueueHandle *command_queue_handle;
#ifdef USE_PROFILING
    command_queue_handle =
    ::dikucl::CommandQueueManager::get_instance().get_command_queue(
                                                                    &err,
                                                                    context_handle,
                                                                    CL_QUEUE_PROFILING_ENABLE);
#else
    command_queue_handle =
    ::dikucl::CommandQueueManager::get_instance().get_command_queue(
                                                                    &err,
                                                                    context_handle);
#endif // USE_PROFILING
    CHECK_CL_ERR(err);
    ::cl::CommandQueue queue = command_queue_handle->command_queue;

    typedef cl_uint KI; const std::string kernel_index_type = "uint";
    typedef cl_float KT; const std::string kernel_real_type = "float";
    typedef cl_float3 KV; const std::string kernel_vector3_type = "float3";

    details::cl::KernelNode<KI, KT, K> *kernel_nodes;
    details::cl::KernelTetrahedron<KI> *kernel_tets;
    details::cl::KernelTetrahedronSurfaceInfo *kernel_tsi;
    details::cl::KernelWorkItem<KI> *kernel_work;
    geometry::ContactsCallback<V> **kernel_callbacks;
    KV *kernel_verts;
    size_t kernel_nodes_size, kernel_tets_size, kernel_verts_size;
    size_t max_bvtt_degree, max_bvtt_height;

    details::cl::gProximityKernelWorkItemGenerator<KI> kernel_work_item_generator;
    details::cl::prepare_input<V, K, T, KI, KT, KV>(
                                                    test_pairs,
                                                    &kernel_nodes, &kernel_nodes_size,
                                                    &kernel_tets, &kernel_tets_size, &kernel_tsi,
                                                    &kernel_verts, &kernel_verts_size,
                                                    &kernel_work_item_generator,
                                                    &kernel_callbacks,
                                                    &max_bvtt_degree, &max_bvtt_height,
                                                    CL_DEVICE_TYPE_GPU, global_mem_cacheline_size,
                                                    (size_t) std::numeric_limits<size_t>::max());

    /* Get Tandem Traversal Kernel and set static parameters */

    ::dikucl::KernelInfo tandem_traversal_kernel_info(
                                                      details::cl::kernels_path + "gproximity/",
                                                      "kdop_cl_gproximity_tandem_traversal.cl",
                                                      "do_tandem_traversal");
    tandem_traversal_kernel_info.define("__K", K);
    tandem_traversal_kernel_info.define("__INDEX_TYPE", kernel_index_type).define("__REAL_TYPE", kernel_real_type);
    tandem_traversal_kernel_info.define("__MAX_BVTT_DEGREE", max_bvtt_degree);
    tandem_traversal_kernel_info.define("__GLOBAL_EXACT_TESTS_CAPACITY", TT_GLOBAL_EXACT_TESTS_CAPACITY);
    tandem_traversal_kernel_info.define("__GLOBAL_WORK_QUEUE_CAPACITY", GLOBAL_WORK_QUEUE_CAPACITY);
    tandem_traversal_kernel_info.define("__IDLE_THRESHOLD", TT_IDLE_THRESHOLD);
    tandem_traversal_kernel_info.define("__LOCAL_WORK_QUEUE_CAPACITY", TT_LOCAL_WORK_QUEUE_CAPACITY);
    tandem_traversal_kernel_info.define("__LOCAL_WORK_QUEUE_INIT_ITEMS", TT_LOCAL_WORK_QUEUE_INIT_ITEMS);
    tandem_traversal_kernel_info.define("__TRAVERSAL_THREADS", TT_THREADS);

    tandem_traversal_kernel_info.include(details::cl::kernels_path).no_signed_zeros(true);
    tandem_traversal_kernel_info.include(details::cl::kernels_path + "gproximity/");
#ifndef NDEBUG
    tandem_traversal_kernel_info.nv_verbose(true);
#endif // NDEBUG
    ::dikucl::KernelHandle *tandem_traversal_kernel_handle = ::dikucl::KernelManager::get_instance().get_kernel(
                                                                                                                tandem_traversal_kernel_info,
                                                                                                                &err,
                                                                                                                context_handle);
    CHECK_CL_ERR(err);
    ::cl::Kernel tandem_traversal_kernel = tandem_traversal_kernel_handle->kernel;

    ::cl::Buffer nodes = ::cl::Buffer(
                                      context, CL_MEM_READ_ONLY,
                                      sizeof (details::cl::KernelNode<KI, KT, K>) * kernel_nodes_size,
                                      NULL, &err);
    CHECK_CL_ERR(err);
    err = queue.enqueueWriteBuffer(
                                   nodes, CL_FALSE,
                                   0, sizeof (details::cl::KernelNode<KI, KT, K>) * kernel_nodes_size,
                                   kernel_nodes);
    CHECK_CL_ERR(err);
    err = tandem_traversal_kernel.setArg(0, nodes);
    CHECK_CL_ERR(err);

    ::cl::Buffer exact_tests = ::cl::Buffer(
                                            context, CL_MEM_READ_WRITE,
                                            TT_GLOBAL_EXACT_TESTS_CAPACITY * sizeof(details::cl::KernelWorkItem<KI>),
                                            NULL, &err);
    CHECK_CL_ERR(err);
    err = tandem_traversal_kernel.setArg(4, exact_tests);
    CHECK_CL_ERR(err);

    /* Get Exact Tests Kernel and set static parameters */

    ::dikucl::KernelInfo exact_tests_kernel_info(
                                                 details::cl::kernels_path,
                                                 "kdop_cl_exact_tests.cl",
                                                 "do_exact_tests");
    exact_tests_kernel_info.define("__INDEX_TYPE", kernel_index_type).define("__REAL_TYPE", kernel_real_type);
    exact_tests_kernel_info.define("__VECTOR3_TYPE", kernel_vector3_type);

    exact_tests_kernel_info.include(details::cl::kernels_path).no_signed_zeros(true);
#ifndef NDEBUG
    exact_tests_kernel_info.nv_verbose(true);
#endif // NDEBUG
    ::dikucl::KernelHandle *exact_tests_kernel_handle = ::dikucl::KernelManager::get_instance().get_kernel(
                                                                                                           exact_tests_kernel_info,
                                                                                                           &err,
                                                                                                           context_handle);
    CHECK_CL_ERR(err);
    ::cl::Kernel exact_tests_kernel = exact_tests_kernel_handle->kernel;

    ::cl::Buffer tetrahedrons = ::cl::Buffer(
                                             context, CL_MEM_READ_ONLY,
                                             sizeof (details::cl::KernelTetrahedron<KI>) * kernel_tets_size,
                                             NULL, &err);
    CHECK_CL_ERR(err);
    err = queue.enqueueWriteBuffer(
                                   tetrahedrons, CL_FALSE,
                                   0, sizeof (details::cl::KernelTetrahedron<KI>) * kernel_tets_size,
                                   kernel_tets);
    CHECK_CL_ERR(err);
    err = exact_tests_kernel.setArg(0, tetrahedrons);
    CHECK_CL_ERR(err);

    ::cl::Buffer tetrahedron_surface_info = ::cl::Buffer(
                                                         context, CL_MEM_READ_ONLY,
                                                         sizeof (details::cl::KernelTetrahedronSurfaceInfo) * kernel_tets_size,
                                                         NULL, &err);
    CHECK_CL_ERR(err);
    err = queue.enqueueWriteBuffer(
                                   tetrahedron_surface_info, CL_FALSE,
                                   0, sizeof (details::cl::KernelTetrahedronSurfaceInfo) * kernel_tets_size,
                                   kernel_tsi);
    CHECK_CL_ERR(err);
    err = exact_tests_kernel.setArg(1, tetrahedron_surface_info);
    CHECK_CL_ERR(err);

    ::cl::Buffer vertices = ::cl::Buffer(
                                         context, CL_MEM_READ_ONLY,
                                         sizeof (KV) * kernel_verts_size,
                                         NULL, &err);
    CHECK_CL_ERR(err);
    err = queue.enqueueWriteBuffer(vertices, CL_FALSE, 0, sizeof (KV) * kernel_verts_size,
                                   kernel_verts);
    CHECK_CL_ERR(err);
    err = exact_tests_kernel.setArg(2, vertices);
    CHECK_CL_ERR(err);

    /* Get Balance Work Kernel and set static parameters */

    ::dikucl::KernelInfo balance_work_kernel_info(
                                                  details::cl::kernels_path + "gproximity/",
                                                  "kdop_cl_gproximity_balance_work.cl",
                                                  "do_balance_work");
    balance_work_kernel_info.define("__INDEX_TYPE", kernel_index_type);
    balance_work_kernel_info.define("__BALANCE_THREADS", BW_THREADS);
    balance_work_kernel_info.define("__GLOBAL_WORK_QUEUE_CAPACITY", GLOBAL_WORK_QUEUE_CAPACITY);
    balance_work_kernel_info.define("__GLOBAL_WORK_QUEUES", TT_WORK_GROUPS);
    balance_work_kernel_info.define("__LOCAL_WORK_QUEUE_INIT_ITEMS", TT_LOCAL_WORK_QUEUE_INIT_ITEMS);
    balance_work_kernel_info.define("__LOCAL_WORK_QUEUE_CAPACITY", TT_LOCAL_WORK_QUEUE_CAPACITY);
    balance_work_kernel_info.define("__TRAVERSAL_THREADS", TT_THREADS);

    balance_work_kernel_info.include(details::cl::kernels_path).no_signed_zeros(true);
    balance_work_kernel_info.include(details::cl::kernels_path + "gproximity/");
#ifndef NDEBUG
    balance_work_kernel_info.nv_verbose(true);
#endif // NDEBUG
    ::dikucl::KernelHandle *balance_work_kernel_handle = ::dikucl::KernelManager::get_instance().get_kernel(
                                                                                                            balance_work_kernel_info,
                                                                                                            &err,
                                                                                                            context_handle);
    CHECK_CL_ERR(err);
    ::cl::Kernel balance_work_kernel = balance_work_kernel_handle->kernel;

    details::cl::KernelContactPoint<KV, KT, KI> *contact_point_container = NULL;
    size_t contact_point_container_size = 0;

    while(!kernel_work_item_generator.empty()) {
      size_t kernel_root_pairs_size = kernel_work_item_generator.generate_kernel_work_items(
                                                                                            TT_WORK_GROUPS, &kernel_work,
                                                                                            global_mem_cacheline_size,
                                                                                            max_bvtt_degree, max_bvtt_height,
                                                                                            TT_WORK_GROUPS,
                                                                                            GLOBAL_WORK_QUEUE_CAPACITY);

      if(kernel_root_pairs_size > 0) {
        ::cl::Buffer work[2];
        work[0] = ::cl::Buffer(
                               context, CL_MEM_READ_WRITE,
                               TT_WORK_GROUPS * GLOBAL_WORK_QUEUE_CAPACITY * sizeof(details::cl::KernelWorkItem<KI>),
                               NULL, &err);
        CHECK_CL_ERR(err);
        err = queue.enqueueWriteBuffer(
                                       work[0], CL_FALSE,
                                       0, TT_WORK_GROUPS * GLOBAL_WORK_QUEUE_CAPACITY * sizeof(details::cl::KernelWorkItem<KI>),
                                       kernel_work);
        CHECK_CL_ERR(err);

        work[1] = ::cl::Buffer(
                               context, CL_MEM_READ_WRITE,
                               TT_WORK_GROUPS * GLOBAL_WORK_QUEUE_CAPACITY * sizeof(details::cl::KernelWorkItem<KI>),
                               NULL, &err);
        CHECK_CL_ERR(err);

        cl_uint *kernel_work_counts =
        new cl_uint[TT_WORK_GROUPS];
        memset(kernel_work_counts, 0,
               TT_WORK_GROUPS * sizeof(cl_uint));
        for(size_t i = 0; i < kernel_root_pairs_size; ++i) {
          kernel_work_counts[i] = 1;
        }
        ::cl::Buffer work_counts = ::cl::Buffer(
                                                context, CL_MEM_READ_WRITE,
                                                TT_WORK_GROUPS * sizeof(cl_uint),
                                                NULL, &err);
        CHECK_CL_ERR(err);
        err = queue.enqueueWriteBuffer(
                                       work_counts, CL_FALSE,
                                       0, TT_WORK_GROUPS * sizeof(cl_uint),
                                       kernel_work_counts);
        CHECK_CL_ERR(err);
        err = tandem_traversal_kernel.setArg(2, work_counts);
        CHECK_CL_ERR(err);

        cl_uint kernel_exact_tests_size = 0;
        ::cl::Buffer exact_tests_size = ::cl::Buffer(
                                                     context, CL_MEM_READ_WRITE,
                                                     sizeof (cl_uint),
                                                     NULL, &err);
        CHECK_CL_ERR(err);
        err = queue.enqueueWriteBuffer(
                                       exact_tests_size, CL_FALSE,
                                       0, sizeof (cl_uint),
                                       &kernel_exact_tests_size);
        CHECK_CL_ERR(err);
        err = tandem_traversal_kernel.setArg(5, exact_tests_size);
        CHECK_CL_ERR(err);

        cl_int kernel_active_splits = 1;
        bool global_work_queue_overflow = false;
        size_t current_work_buffer = 0;
        while (kernel_active_splits > 0 && !global_work_queue_overflow) {
          err = tandem_traversal_kernel.setArg(1, work[current_work_buffer]);
          CHECK_CL_ERR(err);

          cl_uint kernel_idle_count = 0;
          ::cl::Buffer idle_count = ::cl::Buffer(
                                                 context, CL_MEM_READ_WRITE,
                                                 sizeof(cl_uint),
                                                 NULL, &err);
          CHECK_CL_ERR(err);
          err = queue.enqueueWriteBuffer(
                                         idle_count, CL_FALSE,
                                         0, sizeof(cl_uint),
                                         &kernel_idle_count);
          CHECK_CL_ERR(err);
          err = tandem_traversal_kernel.setArg(3, idle_count);
          CHECK_CL_ERR(err);

#ifdef USE_PROFILING
          {
            ::cl::Event tt_profiling;
            err = queue.enqueueNDRangeKernel(
                                             tandem_traversal_kernel,
                                             ::cl::NullRange,
                                             ::cl::NDRange(TT_WORK_GROUPS * TT_THREADS),
                                             ::cl::NDRange(TT_THREADS),
                                             NULL,
                                             &tt_profiling);
            CHECK_CL_ERR(err);
            err = queue.finish();
            CHECK_CL_ERR(err);
            cl_ulong tt_start =
            tt_profiling.getProfilingInfo<CL_PROFILING_COMMAND_START>(&err);
            CHECK_CL_ERR(err);
            cl_ulong tt_end =
            tt_profiling.getProfilingInfo<CL_PROFILING_COMMAND_END>(&err);
            CHECK_CL_ERR(err);
            tandem_traversal_kernel_time += tt_end - tt_start;
          }
#else
          err = queue.enqueueNDRangeKernel(
                                           tandem_traversal_kernel,
                                           ::cl::NullRange,
                                           ::cl::NDRange(TT_WORK_GROUPS * TT_THREADS),
                                           ::cl::NDRange(TT_THREADS));
          CHECK_CL_ERR(err);
#endif // USE_PROFILING
          ++tandem_traversal_invocations;

          err = queue.enqueueReadBuffer(
                                        work_counts, CL_TRUE,
                                        0, TT_WORK_GROUPS * sizeof(cl_uint),
                                        kernel_work_counts);
          CHECK_CL_ERR(err);

          for(size_t i = 0; i < TT_WORK_GROUPS; ++i) {
            if(kernel_work_counts[i] >= GLOBAL_WORK_QUEUE_CAPACITY)
            {
              util::Log logging;

              global_work_queue_overflow = true;

              logging << "Work Queue "
              << i
              << " did overflow ("
              << kernel_work_counts[i]
              << ")"
              << util::Log::newline();
            }
          }

          err = balance_work_kernel.setArg(0, work[current_work_buffer]);
          CHECK_CL_ERR(err);
          err = balance_work_kernel.setArg(1, work[1 - current_work_buffer]);
          CHECK_CL_ERR(err);
          err = balance_work_kernel.setArg(2, work_counts);
          CHECK_CL_ERR(err);

          ::cl::Buffer active_splits = ::cl::Buffer(
                                                    context, CL_MEM_READ_WRITE,
                                                    sizeof(cl_int),
                                                    NULL, &err);
          CHECK_CL_ERR(err);
          err = balance_work_kernel.setArg(3, active_splits);
          CHECK_CL_ERR(err);

          ::cl::Buffer balance_signal = ::cl::Buffer(
                                                     context, CL_MEM_READ_WRITE,
                                                     sizeof(cl_int),
                                                     NULL, &err);
          CHECK_CL_ERR(err);
          err = balance_work_kernel.setArg(4, balance_signal);
          CHECK_CL_ERR(err);

#ifdef USE_PROFILING
          {
            ::cl::Event bw_profiling;
            err = queue.enqueueNDRangeKernel(
                                             balance_work_kernel,
                                             ::cl::NullRange,
                                             ::cl::NDRange(BW_THREADS),
                                             ::cl::NDRange(BW_THREADS),
                                             NULL,
                                             &bw_profiling);
            CHECK_CL_ERR(err);
            err = queue.finish();
            CHECK_CL_ERR(err);
            cl_ulong bw_start =
            bw_profiling.getProfilingInfo<CL_PROFILING_COMMAND_START>(&err);
            CHECK_CL_ERR(err);
            cl_ulong bw_end =
            bw_profiling.getProfilingInfo<CL_PROFILING_COMMAND_END>(&err);
            CHECK_CL_ERR(err);
            balance_work_kernel_time += bw_end - bw_start;
          }
#else
          err = queue.enqueueNDRangeKernel(
                                           balance_work_kernel,
                                           ::cl::NullRange,
                                           ::cl::NDRange(BW_THREADS * BW_WORK_GROUPS),
                                           ::cl::NDRange(BW_THREADS));
          CHECK_CL_ERR(err);
#endif // USE_PROFILING

          err = queue.enqueueReadBuffer(
                                        active_splits, CL_FALSE,
                                        0, sizeof (cl_int),
                                        &kernel_active_splits);
          CHECK_CL_ERR(err);

          cl_int kernel_balance_signal = 0;
          err = queue.enqueueReadBuffer(
                                        balance_signal, CL_TRUE,
                                        0, sizeof (cl_int),
                                        &kernel_balance_signal);
          CHECK_CL_ERR(err);

          if(kernel_balance_signal == 1) {
            current_work_buffer = 1 - current_work_buffer;
          }
        }

        delete[] kernel_work_counts;

        err = queue.enqueueReadBuffer(
                                      exact_tests_size, CL_TRUE,
                                      0, sizeof (cl_uint),
                                      &kernel_exact_tests_size);
        CHECK_CL_ERR(err);

        if (kernel_exact_tests_size > 0) {
          err = exact_tests_kernel.setArg(3, exact_tests);
          CHECK_CL_ERR(err);
          err = exact_tests_kernel.setArg(4, kernel_exact_tests_size);
          CHECK_CL_ERR(err);

          ::cl::Buffer contact_points = ::cl::Buffer(
                                                     context, CL_MEM_WRITE_ONLY,
                                                     sizeof (details::cl::KernelContactPoint<KV, KT, KI>) * kernel_exact_tests_size * 20,
                                                     NULL, &err);
          CHECK_CL_ERR(err);
          err = exact_tests_kernel.setArg(5, contact_points);
          CHECK_CL_ERR(err);

          cl_uint kernel_contact_points_size = 0;
          ::cl::Buffer contact_points_size = ::cl::Buffer(
                                                          context, CL_MEM_READ_WRITE,
                                                          sizeof (cl_uint),
                                                          NULL, &err);
          CHECK_CL_ERR(err);
          err = queue.enqueueWriteBuffer(
                                         contact_points_size, CL_FALSE,
                                         0, sizeof (cl_uint),
                                         &kernel_contact_points_size);
          CHECK_CL_ERR(err);
          err = exact_tests_kernel.setArg(6, contact_points_size);
          CHECK_CL_ERR(err);

#ifdef USE_PROFILING
          {
            ::cl::Event et_profiling;
            err = queue.enqueueNDRangeKernel(
                                             exact_tests_kernel,
                                             ::cl::NullRange,
                                             ::cl::NDRange(ET_THREADS * ET_WORK_GROUPS),
                                             ::cl::NDRange(ET_THREADS),
                                             NULL,
                                             &et_profiling);
            CHECK_CL_ERR(err);
            err = queue.finish();
            CHECK_CL_ERR(err);
            cl_ulong et_start = et_profiling.getProfilingInfo<CL_PROFILING_COMMAND_START>(&err);
            CHECK_CL_ERR(err);
            cl_ulong et_end = et_profiling.getProfilingInfo<CL_PROFILING_COMMAND_END>(&err);
            CHECK_CL_ERR(err);
            exact_tests_kernel_time += et_end - et_start;
          }
#else
          err = queue.enqueueNDRangeKernel(
                                           exact_tests_kernel,
                                           ::cl::NullRange,
                                           ::cl::NDRange(ET_THREADS * ET_WORK_GROUPS),
                                           ::cl::NDRange(ET_THREADS));
          CHECK_CL_ERR(err);
#endif // USE_PROFILING
          ++exact_test_invocations;

          err = queue.enqueueReadBuffer(
                                        contact_points_size, CL_TRUE,
                                        0, sizeof (cl_uint),
                                        &kernel_contact_points_size);
          CHECK_CL_ERR(err);

          if (kernel_contact_points_size > 0) {
            contact_point_container = (details::cl::KernelContactPoint<KV, KT, KI>*) realloc(
                                                                                             contact_point_container,
                                                                                             sizeof (details::cl::KernelContactPoint<KV, KT, KI>) *
                                                                                             (contact_point_container_size + (size_t) kernel_contact_points_size));

            err = queue.enqueueReadBuffer(
                                          contact_points, CL_TRUE,
                                          0, sizeof (details::cl::KernelContactPoint<KV, KT, KI>) * (size_t) kernel_contact_points_size,
                                          contact_point_container + (size_t) contact_point_container_size);
            CHECK_CL_ERR(err);
            contact_point_container_size += (size_t) kernel_contact_points_size;                        
          }
          
        }
        
      }
      
      kernel_work_item_generator.cleanup_generated_work_items(
                                                              kernel_work);
    }
    
    err = queue.finish();
    CHECK_CL_ERR(err);
    
    if(contact_point_container != NULL) {
      qsort(  contact_point_container,
            contact_point_container_size,
            sizeof(details::cl::KernelContactPoint<KV, KT, KI>),
            details::cl::compare_kernel_contact_points<KV, KT, KI>);
      for(size_t i = 0; i < contact_point_container_size; ++i) {
        details::cl::KernelContactPoint<KV, KT, KI> contact_point = contact_point_container[i];
        V p = V::make((T) contact_point.p.s[0], (T) contact_point.p.s[1], (T) contact_point.p.s[2]);
        V n = V::make((T) contact_point.n.s[0], (T) contact_point.n.s[1], (T) contact_point.n.s[2]);
        T d = (T) contact_point.d;
        kernel_callbacks[(size_t) contact_point.tp]->operator()(p, n, d);
      }
      free(contact_point_container);
    }
    
    details::cl::cleanup<V, K, T, KI, KT, KV>(
                                              kernel_nodes,
                                              kernel_tets,
                                              kernel_tsi,
                                              kernel_verts,
                                              kernel_callbacks,
                                              CL_DEVICE_TYPE_GPU);
    
    details::cl::record_kernel_times(  tandem_traversal_kernel_time + balance_work_kernel_time
                                     , exact_tests_kernel_time);
    details::cl::record_kernel_invocations(  tandem_traversal_invocations
                                           , exact_test_invocations);
  }
  
} // namespace kdop

// KDOP_CL_GPROXIMITY_TANDEM_TRAVERSAL_H
#endif
