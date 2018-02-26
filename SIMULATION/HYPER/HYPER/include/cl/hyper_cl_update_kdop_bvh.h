#ifndef HYPER_CL_UPDATE_KDOP_BVH_H
#define HYPER_CL_UPDATE_KDOP_BVH_H

#include <cl/hyper_cl_kernels_path.h>
#include <cl/kdop_cl_refit_tree.h>
#include <cl/kdop_cl_tree.h>

#include <dikucl.hpp>
#include <dikucl_command_queue_manager.hpp>
#include <dikucl_context_manager.hpp>
#include <dikucl_kernel_manager.hpp>
#include <dikucl_util.h>

#include <hyper_body.h>
#include <hyper_engine.h>

#include <kdop_tree.h>

#include <mesh_array.h>

#include <util_profiling.h>

#include <cassert>
#include <string>
#include <vector>

namespace hyper
{
    
  template<typename MT>
  inline void update_kdop_bvh(
                              Engine<MT> & engine
                              , kdop::dikucl const & /* tag */
                              , size_t open_cl_platform = 0
                              , size_t open_cl_device = 0
                              )
  {
    typedef typename MT::real_type     T;
    typedef typename MT::vector3_type  V;
    typedef tiny::ValueTraits<T>      VT;
    const size_t K = 8;
    
    typedef typename std::vector< Body< MT > >::iterator object_iterator;
    
    // Kernel types
    typedef cl_float  KT;
    typedef cl_float3 KV;
    typedef cl_float4 KQ;

    cl_ulong refit_tree_kernel_time       = 0;
    cl_ulong refit_tree_opencl_read_time  = 0;
    cl_ulong refit_tree_opencl_write_time = 0;
    
    cl_int err = CL_SUCCESS;

    START_TIMER("update_kdop_opencl_init");
    
    // Specified platform
    ::dikucl::PlatformHandle *platform_handle =
            ::dikucl::PlatformManager::get_instance().get_platform(  &err
                                                                   , open_cl_platform);
    CHECK_CL_ERR(err);
    
    // Specified device
    ::dikucl::DeviceHandle *device_handle =
            ::dikucl::DeviceManager::get_instance().get_device(  &err
                                                             , platform_handle
                                                             , open_cl_device
                                                             );
    CHECK_CL_ERR(err);
    cl_device_type device_type = device_handle->device.getInfo<CL_DEVICE_TYPE>(&err);
    CHECK_CL_ERR(err);
    
    // Cacheline size is needed when we're on the CPU in order to ensure proper
    // memory alignment.
    size_t global_mem_cacheline_size =
            (size_t) device_handle->device.getInfo<CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE>(&err);
    CHECK_CL_ERR(err);
    
    // Default context
    ::dikucl::ContextHandle *context_handle =
            ::dikucl::ContextManager::get_instance().get_context(&err, device_handle);
    CHECK_CL_ERR(err);
    ::cl::Context context = context_handle->context;
    
    // Get command queue for that context
    ::dikucl::CommandQueueHandle *command_queue_handle;
#ifdef USE_PROFILING
    command_queue_handle =
            ::dikucl::CommandQueueManager::get_instance().get_command_queue(  &err
                                                                          , context_handle
                                                                          , CL_QUEUE_PROFILING_ENABLE
                                                                          );
#else
    command_queue_handle =
            ::dikucl::CommandQueueManager::get_instance().get_command_queue(  &err
                                                                          , context_handle
                                                                          );
#endif // USE_PROFILING            
    CHECK_CL_ERR(err);
    ::cl::CommandQueue queue = command_queue_handle->command_queue;
    STOP_TIMER("update_kdop_opencl_init");

    START_TIMER("refit_tree_opencl_prepare_input");
    
    object_iterator current   = engine.body_begin();
    const object_iterator end = engine.body_end();
    
    /* Gather information for all objects */
    
    // Lookup table for each chunk of each level of each object:
    // i.e. when flattened out (= the levels and chunks are concatenated into one
    // vector, starting with the topmost level/leftmost chunk), what would the
    // offsets for each level and chunk be. [0][0][0] will tell the offset
    // of the first chunk of the first level of the first object (which will be 0
    // for each separate object). [1][2][3] tells the offset of the 4th chunk
    // in the 3rd level of the 2nd object and so forth.
    // This information is needed to lay out memory properly on the OpenCL
    // device.
    size_t ***flattened_object_node_offsets = new size_t**[std::distance(current, end)];
    
    // Tracks the overall number of objects
    size_t total_objects = 0;
    
    // Tracks the overall number of vertices, tetrahedrons and nodes.
    // Again, these numbers are needed to lay out memory properly on the OpenCL
    // device.
    size_t total_vertices = 0, total_tetrahedrons = 0, total_nodes = 0;
    for (; current != end; ++current, ++total_objects)
    {
      if(current->empty()) {
        continue;
      }
      
      total_vertices     += current->m_mesh.vertex_size();
      total_tetrahedrons += current->m_mesh.tetrahedron_size();
      
      // Offsets are not global but per-object. I.e. the first chunk in the first
      // level of each object will always have offset 0.
      size_t total_nodes_per_object = 0;
      
      flattened_object_node_offsets[total_objects] = new size_t*[current->m_tree.number_of_levels()];
      for(size_t h = 0; h < current->m_tree.number_of_levels(); ++h)
      {
        flattened_object_node_offsets[total_objects][h] = new size_t[current->m_tree.super_chunks(h).size()];
        for(size_t c = 0; c < current->m_tree.super_chunks(h).size(); ++c)
        {
          // The offset of the current chunk in the current level is the number
          // of nodes encountered so far.
          flattened_object_node_offsets[total_objects][h][c] = total_nodes_per_object;
          kdop::SubTree<T, K> chunk = current->m_tree.super_chunks(h)[c];
          for(size_t n = 0; n < chunk.m_nodes.size(); ++n)
          {
            if(! chunk.m_nodes[n].is_undefined()) {
              ++total_nodes;
              ++total_nodes_per_object;
            }
          }
        }
      }
    }
    
    // Make room for kernel data.
    
    // Holds the starting offset of vertices for each object. This is needed
    // because all vertices for all objects will reside in one array.
    cl_uint *vertex_offsets                                     = NULL;
    
    // Repeat for tetrahedrons and BVH nodes.
    cl_uint *tetrahedron_offsets                                = NULL;
    cl_uint *node_offsets                                       = NULL;
    
    // Holds the offsets of the last chunk level (i.e. the first non-
    // superchunk of the object). This is needed later in the kernel to
    // identify whether an encountered leaf is a leaf of a super-chunk pointing
    // to another chunk, or a leaf of a branch pointing to a tetrahedron.
    cl_uint *last_level_offsets                                 = NULL;
    
    // The actual data.
    KV *vertices                                                = NULL;
    kdop::details::cl::KernelTetrahedron<cl_uint> *tetrahedrons = NULL;
    kdop::details::cl::KernelNode<cl_uint, KT, K> *nodes        = NULL;
    
    if (device_type & CL_DEVICE_TYPE_CPU)
    {
      // On the CPU we need to align memory properly so the device can use host
      // pointers without overhead.
      vertex_offsets      = (cl_uint*) ALIGNED_ALLOC(4096, total_objects * sizeof(cl_uint));
      tetrahedron_offsets = (cl_uint*) ALIGNED_ALLOC(4096, total_objects * sizeof(cl_uint));
      node_offsets        = (cl_uint*) ALIGNED_ALLOC(4096, total_objects * sizeof(cl_uint));
      last_level_offsets  = (cl_uint*) ALIGNED_ALLOC(4096, total_objects * sizeof(cl_uint));
      vertices            = (KV*) ALIGNED_ALLOC(4096, total_vertices * sizeof(KV));
      tetrahedrons        = (kdop::details::cl::KernelTetrahedron<cl_uint>*) ALIGNED_ALLOC(
              4096, total_tetrahedrons * sizeof(kdop::details::cl::KernelTetrahedron<cl_uint>)
              );
      nodes               = (kdop::details::cl::KernelNode<cl_uint, KT, K>*) ALIGNED_ALLOC(
              4096, total_nodes * sizeof(kdop::details::cl::KernelNode<cl_uint, KT, K>)
              );
    } else {
      vertex_offsets      = new cl_uint[total_objects];
      tetrahedron_offsets = new cl_uint[total_objects];
      node_offsets        = new cl_uint[total_objects];
      last_level_offsets  = new cl_uint[total_objects];
      vertices            = new KV[total_vertices];
      tetrahedrons        = new kdop::details::cl::KernelTetrahedron<cl_uint>[total_tetrahedrons];
      nodes               = new kdop::details::cl::KernelNode<cl_uint, KT, K>[total_nodes];
    }
    
    /* Initialize kernel data. */
    
    // Reset counters because we go through all objects again now.
    current = engine.body_begin();
    total_objects = 0;
    total_vertices = 0;
    total_tetrahedrons = 0;
    total_nodes = 0;
    for (; current != end; ++current, ++total_objects)
    {
      if(current->empty()) {
        continue;
      }
            
      size_t N = current->m_mesh.vertex_size();
      for(size_t n = 0u; n < N; ++n)
      {
        mesh_array::Vertex const & v = current->m_mesh.vertex(n);
        
        // Write vertex to appropriate position. Vertex offset for the first
        // object is 0, the other ones are the cumulated counts of vertices of
        // the objects before them (see directly below).
        const size_t vertex_offset = (total_objects == 0 ? 0 : vertex_offsets[total_objects - 1]) + v.idx();
        vertices[vertex_offset].s[0] = current->m_X(v);
        vertices[vertex_offset].s[1] = current->m_Y(v);
        vertices[vertex_offset].s[2] = current->m_Z(v);
      }
      
      // Remember the number of vertices per object (cumulative).
      total_vertices += N;
      vertex_offsets[total_objects] = total_vertices;
      
      // Repeat for tetrahedrons.
      N = current->m_mesh.tetrahedron_size();
      for(size_t n = 0u; n < N; ++n)
      {
        mesh_array::Tetrahedron const & t = current->m_mesh.tetrahedron(n);
          
        const size_t tetrahedron_offset = (total_objects == 0 ? 0 : tetrahedron_offsets[total_objects - 1]) + t.idx();
        tetrahedrons[tetrahedron_offset].vertex_idx[0] = t.i();
        tetrahedrons[tetrahedron_offset].vertex_idx[1] = t.j();
        tetrahedrons[tetrahedron_offset].vertex_idx[2] = t.k();
        tetrahedrons[tetrahedron_offset].vertex_idx[3] = t.m();
      }
      
      total_tetrahedrons += N;
      tetrahedron_offsets[total_objects] = total_tetrahedrons;
      
      // Repeat for BVH nodes.
      size_t total_nodes_per_object = 0;
      for(size_t h = 0; h < current->m_tree.number_of_levels(); ++h)
      {                
        for(size_t c = 0; c < current->m_tree.super_chunks(h).size(); ++c)
        {                    
          kdop::SubTree<T, K> chunk = current->m_tree.super_chunks(h)[c];
          
          for(size_t n = 0; n < chunk.m_nodes.size(); ++n)
          {
            kdop::Node<T, K> const & node = chunk.m_nodes[n];
            
            if(!node.is_undefined()) {
              const size_t node_offset = (total_objects == 0 ? 0 : node_offsets[total_objects - 1]) + total_nodes_per_object;
                            
              if(node.is_leaf() && h == current->m_tree.number_of_levels() - 1) {
                // Index into tetrahedrons, no need to offset because the kernel
                // will use the tetrahedron_offsets.
                nodes[node_offset].start = node.m_start;
                nodes[node_offset].end   = node.m_end;
              } else if(node.is_leaf()) {
                // Index into next level, so offset by the offset of the next
                // level.
                nodes[node_offset].start = flattened_object_node_offsets[total_objects][h + 1][node.m_start];
                nodes[node_offset].end   = flattened_object_node_offsets[total_objects][h + 1][node.m_end];
              } else {
                // Index into this chunk, so offset by the offset of this chunk
                // in this level.
                nodes[node_offset].start = flattened_object_node_offsets[total_objects][h][c] + node.m_start;
                nodes[node_offset].end   = flattened_object_node_offsets[total_objects][h][c] + node.m_end;
              }
              
              for(size_t k = 0; k < K / 2; ++k)
              {
                nodes[node_offset].slabs[k].lower = VT::highest();
                nodes[node_offset].slabs[k].upper = VT::lowest();
              }

              ++total_nodes;
              ++total_nodes_per_object;
            }
          }
        }
      }
      
      node_offsets[total_objects] = total_nodes;
      
      // Remember the index of the first chunk in the last level. This is so
      // the kernel can identify what kind of leaf it deals with (i.e. a leaf
      // referencing another chunk or a tetrahedron).
      last_level_offsets[total_objects] = flattened_object_node_offsets[total_objects][current->m_tree.number_of_levels() - 1][0];
    }
    
    // Intermediate cleanup.
    current = engine.body_begin();
    total_objects = 0;
    for(; current != end; ++current, ++total_objects) {
        if(current->empty()) {
            continue;
        }
        for(size_t h = 0; h < current->m_tree.number_of_levels(); ++h) {
            delete[] flattened_object_node_offsets[total_objects][h];
        }
        delete[] flattened_object_node_offsets[total_objects];
    }
    delete[] flattened_object_node_offsets;

    STOP_TIMER("refit_tree_opencl_prepare_input");
    START_TIMER("update_kdop_opencl_host_time");
    
    /* Copy over data and set kernel arguments. */
    
    ::cl::Buffer vertices_buffer;
#ifdef USE_PROFILING
    ::cl::Event vertices_buffer_event;
#endif // USE_PROFILING
    if (device_type & CL_DEVICE_TYPE_CPU) {
      size_t size = total_vertices * sizeof(KV);
      size += size % global_mem_cacheline_size;
      vertices_buffer = ::cl::Buffer(  context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR
                                     , size
                                     , vertices, &err
                                     );
      CHECK_CL_ERR(err);
    } else {
      vertices_buffer = ::cl::Buffer(  context, CL_MEM_READ_WRITE
                                     , total_vertices * sizeof(KV)
                                     , NULL, &err
                                     );
      CHECK_CL_ERR(err);
      err = queue.enqueueWriteBuffer(  vertices_buffer, CL_FALSE
                                     , 0, total_vertices * sizeof(KV)
                                     , vertices
                                     , NULL,
#ifdef USE_PROFILING
                                     &vertices_buffer_event);
#else // USE_PROFILING
                                     NULL);
#endif // USE_PROFILING
      CHECK_CL_ERR(err);
    }
    
    ::cl::Buffer vertex_offsets_buffer;
#ifdef USE_PROFILING
    ::cl::Event vertex_offsets_buffer_event;
#endif // USE_PROFILING
    if (device_type & CL_DEVICE_TYPE_CPU) {
      size_t size = total_objects * sizeof(cl_uint);
      size += size % global_mem_cacheline_size;
      vertex_offsets_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR
                                           , size
                                           , vertex_offsets, &err
                                      );
      CHECK_CL_ERR(err);
    } else {
      vertex_offsets_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY
                                           , total_objects * sizeof(cl_uint)
                                           , NULL, &err
                                      );
      CHECK_CL_ERR(err);
      err = queue.enqueueWriteBuffer(  vertex_offsets_buffer, CL_FALSE
                                     , 0, total_objects * sizeof(cl_uint)
                                     , vertex_offsets
                                     , NULL,
#ifdef USE_PROFILING
                                     &vertex_offsets_buffer_event);
#else // USE_PROFILING
                                     NULL);
#endif // USE_PROFILING
      CHECK_CL_ERR(err);
    }
    
    /* Prepare data and buffers for the refitting kernel. */
    
    ::cl::Buffer tetrahedrons_buffer;
#ifdef USE_PROFILING
    ::cl::Event tetrahedrons_buffer_event;
#endif // USE_PROFILING
    if (device_type & CL_DEVICE_TYPE_CPU) {
      size_t size = total_tetrahedrons * sizeof(kdop::details::cl::KernelTetrahedron<cl_uint>);
      size += size % global_mem_cacheline_size;
      tetrahedrons_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR
                                         , size
                                         , tetrahedrons, &err
                                         );
      CHECK_CL_ERR(err);
    } else {
      tetrahedrons_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY
                                         , total_tetrahedrons * sizeof(kdop::details::cl::KernelTetrahedron<cl_uint>)
                                         , NULL, &err
                                         );
      CHECK_CL_ERR(err);
      err = queue.enqueueWriteBuffer(  tetrahedrons_buffer, CL_FALSE
                                     , 0, total_tetrahedrons * sizeof(kdop::details::cl::KernelTetrahedron<cl_uint>)
                                     , tetrahedrons
                                     , NULL,
#ifdef USE_PROFILING
                                     &tetrahedrons_buffer_event);
#else // USE_PROFILING
                                     NULL);
#endif // USE_PROFILING
      CHECK_CL_ERR(err);
    }
    
    ::cl::Buffer tetrahedron_offsets_buffer;
#ifdef USE_PROFILING
    ::cl::Event tetrahedron_offsets_buffer_event;
#endif // USE_PROFILING
    if (device_type & CL_DEVICE_TYPE_CPU) {
      size_t size = total_objects * sizeof(cl_uint);
      size += size % global_mem_cacheline_size;
      tetrahedron_offsets_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR
                                                , size
                                                , tetrahedron_offsets, &err
                                                );
      CHECK_CL_ERR(err);
    } else {
      tetrahedron_offsets_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY
                                                , total_objects * sizeof(cl_uint)
                                                , NULL, &err
                                                );
      CHECK_CL_ERR(err);
      err = queue.enqueueWriteBuffer(  tetrahedron_offsets_buffer, CL_FALSE
                                     , 0, total_objects * sizeof(cl_uint)
                                     , tetrahedron_offsets
                                     , NULL,
#ifdef USE_PROFILING
                                     &tetrahedron_offsets_buffer_event);
#else // USE_PROFILING
                                     NULL);
#endif // USE_PROFILING
      CHECK_CL_ERR(err);
    }
    
    // The nodes_buffer is for input and output as it will receive the updated
    // BVH nodes.
    ::cl::Buffer nodes_buffer;
#ifdef USE_PROFILING
    ::cl::Event write_nodes_buffer_event;
#endif // USE_PROFILING
    if (device_type & CL_DEVICE_TYPE_CPU) {
      size_t size = total_nodes * sizeof(kdop::details::cl::KernelNode<cl_uint, KT, K>);
      size += size % global_mem_cacheline_size;
      nodes_buffer = ::cl::Buffer(  context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR
                                  , size
                                  , nodes, &err
                                  );
      CHECK_CL_ERR(err);
    } else {
      nodes_buffer = ::cl::Buffer(  context, CL_MEM_READ_WRITE
                                         , total_nodes * sizeof(kdop::details::cl::KernelNode<cl_uint, KT, K>)
                                         , NULL, &err
                                         );
      CHECK_CL_ERR(err);
      err = queue.enqueueWriteBuffer(  nodes_buffer, CL_FALSE
                                     , 0, total_nodes * sizeof(kdop::details::cl::KernelNode<cl_uint, KT, K>)
                                     , nodes
                                     , NULL,
#ifdef USE_PROFILING
                                     &write_nodes_buffer_event);
#else // USE_PROFILING
                                     NULL);
#endif // USE_PROFILING
      CHECK_CL_ERR(err);
    }
    
    ::cl::Buffer node_offsets_buffer;
    if (device_type & CL_DEVICE_TYPE_CPU) {
      size_t size = total_objects * sizeof(cl_uint);
      size += size % global_mem_cacheline_size;
      node_offsets_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR
                                         , size
                                         , node_offsets, &err
                                         );
      CHECK_CL_ERR(err);
    } else {
      node_offsets_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY
                                         , total_objects * sizeof(cl_uint)
                                         , NULL, &err
                                         );
      CHECK_CL_ERR(err);
      err = queue.enqueueWriteBuffer(  node_offsets_buffer, CL_FALSE
                                     , 0, total_objects * sizeof(cl_uint)
                                     , node_offsets
                                     );
      CHECK_CL_ERR(err);
    }
    
    ::cl::Buffer last_level_offsets_buffer;
#ifdef USE_PROFILING
    ::cl::Event last_level_offsets_buffer_event;
#endif // USE_PROFILING
    if (device_type & CL_DEVICE_TYPE_CPU) {
      size_t size = total_objects * sizeof(cl_uint);
      size += size % global_mem_cacheline_size;
      last_level_offsets_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR
                                               , size
                                               , last_level_offsets, &err
                                               );
      CHECK_CL_ERR(err);
    } else {
      last_level_offsets_buffer = ::cl::Buffer(  context, CL_MEM_READ_ONLY
                                               , total_objects * sizeof(cl_uint)
                                               , NULL, &err
                                               );
      CHECK_CL_ERR(err);
      err = queue.enqueueWriteBuffer(  last_level_offsets_buffer, CL_FALSE
                                     , 0, total_objects * sizeof(cl_uint)
                                     , last_level_offsets
                                     , NULL,
#ifdef USE_PROFILING
                                     &last_level_offsets_buffer_event);
#else // USE_PROFILING
                                     NULL);
#endif // USE_PROFILING
      CHECK_CL_ERR(err);
    }
    
    // Invoke internals of the KDOP module.
    refit_tree_kernel_time +=
        kdop::details::cl::refit_tree<V, K, T>(  context_handle, device_handle, &queue
                                               , &tetrahedrons_buffer, &tetrahedron_offsets_buffer
                                               , &vertices_buffer, &vertex_offsets_buffer
                                               , &nodes_buffer, &node_offsets_buffer, &last_level_offsets_buffer
                                               , total_objects
                                               );
    
    // Queue reading back of refitted BVH nodes.
#ifdef USE_PROFILING
    ::cl::Event read_nodes_buffer_event;
#endif // USE_PROFILING
    if(device_type & CL_DEVICE_TYPE_CPU) {
      kdop::details::cl::KernelNode<cl_uint, KT, K> *mapped_nodes =
          (kdop::details::cl::KernelNode<cl_uint, KT, K>*) queue.enqueueMapBuffer(
                nodes_buffer, CL_FALSE, CL_MAP_READ
              , 0, total_nodes * sizeof(kdop::details::cl::KernelNode<cl_uint, KT, K>)
              , NULL, NULL, &err
              );
      CHECK_CL_ERR(err);
      std::memmove(nodes, mapped_nodes, total_nodes * sizeof(kdop::details::cl::KernelNode<cl_uint, KT, K>));
      err = queue.enqueueUnmapMemObject(nodes_buffer, mapped_nodes);
      CHECK_CL_ERR(err);
    } else {
      err = queue.enqueueReadBuffer(
              nodes_buffer, CL_FALSE
              , 0, total_nodes * sizeof(kdop::details::cl::KernelNode<cl_uint, KT, K>)
              , nodes
              , NULL,
#ifdef USE_PROFILING
              &read_nodes_buffer_event);
#else // USE_PROFILING
              NULL);
#endif // USE_PROFILING
      CHECK_CL_ERR(err);
    }
    
    // Wait for all pending operations to finish.
    err = queue.finish();
    CHECK_CL_ERR(err);

#ifdef USE_PROFILING
    {
      if((device_type & CL_DEVICE_TYPE_CPU) == 0) {
        // on the CPU we use host memory, no transactions to profile
        ::cl::Event refit_tree_write_events[] = { tetrahedrons_buffer_event, tetrahedron_offsets_buffer_event,
          write_nodes_buffer_event, last_level_offsets_buffer_event
        };
        for (short i = 0; i < 4; ++i) {
          cl_ulong start = refit_tree_write_events[i].getProfilingInfo<CL_PROFILING_COMMAND_START>(&err);
          CHECK_CL_ERR(err);
          cl_ulong end = refit_tree_write_events[i].getProfilingInfo<CL_PROFILING_COMMAND_END>(&err);
          CHECK_CL_ERR(err);
          refit_tree_opencl_write_time += end - start;
        }

        cl_ulong start = read_nodes_buffer_event.getProfilingInfo<CL_PROFILING_COMMAND_START>(&err);
        CHECK_CL_ERR(err);
        cl_ulong end = read_nodes_buffer_event.getProfilingInfo<CL_PROFILING_COMMAND_END>(&err);
        CHECK_CL_ERR(err);
        refit_tree_opencl_read_time += end - start;
      }
    }
#endif // USE_PROFILING
    STOP_TIMER("update_kdop_opencl_host_time");

    START_TIMER("refit_tree_opencl_prepare_output");
    
    // Run through all objects again and update them accordingly.
    current = engine.body_begin();
    total_objects = 0;
    for (; current != end; ++current, ++total_objects)
    {      
      if(current->empty()) {
          continue;
      }
      
      // Use the same calculations as above for figuring out all offsets.
      size_t total_nodes_per_object = 0;
      for(size_t h = 0; h < current->m_tree.number_of_levels(); ++h)
      {
        for(size_t c = 0; c < current->m_tree.super_chunks(h).size(); ++c)
        {
          kdop::SubTree<T, K> & chunk = current->m_tree.super_chunks(h)[c];
          for(size_t n = 0; n < chunk.m_nodes.size(); ++n)
          {
            if(! chunk.m_nodes[n].is_undefined())
            {
              kdop::Node<T, K> & node = chunk.m_nodes[n];
              const size_t node_offset = (total_objects == 0 ? 0 : node_offsets[total_objects - 1]) + total_nodes_per_object;
              
              for(size_t k = 0; k < K / 2; ++k)
              {
                  node.m_volume(k).lower() = (T) nodes[node_offset].slabs[k].lower;
                  node.m_volume(k).upper() = (T) nodes[node_offset].slabs[k].upper;
              }
              
              ++total_nodes_per_object;
            }
          }
        }
      }
    }
    
    // Cleanup.
    if (device_type & CL_DEVICE_TYPE_CPU)
    {
      free(vertex_offsets);
      free(tetrahedron_offsets);
      free(node_offsets);
      free(last_level_offsets);
      free(vertices);
      free(tetrahedrons);
      free(nodes);
    } else {
      delete[] vertex_offsets;
      delete[] tetrahedron_offsets;
      delete[] node_offsets;
      delete[] last_level_offsets;
      delete[] vertices;
      delete[] tetrahedrons;
      delete[] nodes;
    }

    STOP_TIMER("refit_tree_opencl_prepare_output");
    
    RECORD_TIME("refit_tree", (double) refit_tree_kernel_time / 1000000.0);
    RECORD_TIME("refit_tree_opencl_read_time", refit_tree_opencl_read_time / 1000000.0);
    RECORD_TIME("refit_tree_opencl_write_time", refit_tree_opencl_write_time / 1000000.0);
    RECORD("refit_tree_invocations", 1);
  }
    
} // end of namespace hyper

#endif // HYPER_CL_UPDATE_KDOP_BVH_H