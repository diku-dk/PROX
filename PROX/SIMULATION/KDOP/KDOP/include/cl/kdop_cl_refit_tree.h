#ifndef KDOP_CL_REFIT_TREE_H
#define KDOP_CL_REFIT_TREE_H

#include <cassert>

#include <dikucl.hpp>
#include <dikucl_context_manager.hpp>
#include <dikucl_device_manager.hpp>
#include <dikucl_kernel_manager.hpp>
#include <dikucl_util.h>

#include <kdop_tree.h>

#include <mesh_array_t4mesh.h>
#include <mesh_array_vertex_attribute.h>

#include "kdop_cl_kernels_path.h"

namespace kdop
{
    namespace details
    {
        
        namespace cl
        {
    
            /**
             * Starts the refitting BVH kernel, assuming all data is scheduled
             * to be copied to the device in queue. Does not finish the queue.
             * 
             * @param context_handle
             * @param device_handle
             * @param queue
             * @param tetrahedrons
             * @param tetrahedron_offsets
             * @param vertices
             * @param vertex_offsets
             * @param nodes
             * @param node_offsets
             * @param last_level_offsets
             * @param total_objects
             */
            template<typename V, size_t K, typename T>
            inline cl_ulong refit_tree(  ::dikucl::ContextHandle *context_handle
                                       , ::dikucl::DeviceHandle *device_handle
                                       , ::cl::CommandQueue *queue
                                       , ::cl::Buffer *tetrahedrons
                                       , ::cl::Buffer *tetrahedron_offsets
                                       , ::cl::Buffer *vertices
                                       , ::cl::Buffer *vertex_offsets
                                       , ::cl::Buffer *nodes /* out */
                                       , ::cl::Buffer *node_offsets
                                       , ::cl::Buffer *last_level_offsets
                                       , size_t total_objects
                                   )
            {
                cl_int err = CL_SUCCESS;
                
                // Get the BVH refitting kernel.
                ::dikucl::KernelInfo refit_tree_kernel_info(  details::cl::kernels_path
                                                          , "kdop_cl_refit_tree.cl"
                                                          , "do_refit_tree");
                refit_tree_kernel_info.include(details::cl::kernels_path).no_signed_zeros(true);
                refit_tree_kernel_info.define("__K", K);
                refit_tree_kernel_info.define("__INDEX_TYPE", "uint");
#ifndef NDEBUG
                refit_tree_kernel_info.nv_verbose(true).debug_intel(true);
#endif // NDEBUG
                ::dikucl::KernelHandle *refit_tree_kernel_handle =
                        ::dikucl::KernelManager::get_instance().get_kernel(  refit_tree_kernel_info
                                                                         , &err, context_handle);
                CHECK_CL_ERR(err);
                ::cl::Kernel refit_tree_kernel = refit_tree_kernel_handle->kernel;
                
                // Simply set all arguments.
                err = refit_tree_kernel.setArg(0, *vertices);
                CHECK_CL_ERR(err);
                err = refit_tree_kernel.setArg(1, *vertex_offsets);
                CHECK_CL_ERR(err);
                err = refit_tree_kernel.setArg(2, *tetrahedrons);
                CHECK_CL_ERR(err);
                err = refit_tree_kernel.setArg(3, *tetrahedron_offsets);
                CHECK_CL_ERR(err);
                err = refit_tree_kernel.setArg(4, *nodes);
                CHECK_CL_ERR(err);
                err = refit_tree_kernel.setArg(5, *node_offsets);
                CHECK_CL_ERR(err);
                err = refit_tree_kernel.setArg(6, *last_level_offsets);
                CHECK_CL_ERR(err);
                
                // Use system suggested local work size.
                const size_t rt_local_work_size =
                        refit_tree_kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(  device_handle->device
                                                                                      , &err);
                CHECK_CL_ERR(err);
                
                // Start the kernel.
                cl_ulong duration = 0;
#ifdef USE_PROFILING
                {
                    ::cl::Event rt_profiling;
                    err = queue->enqueueNDRangeKernel(  refit_tree_kernel
                                                     , ::cl::NullRange
                                                     , ::cl::NDRange(total_objects * rt_local_work_size)
                                                     , ::cl::NDRange(rt_local_work_size)
                                                     , NULL
                                                     , &rt_profiling
                                                     );
                    CHECK_CL_ERR(err);
                    err = queue->finish();
                    CHECK_CL_ERR(err);
                    cl_ulong start = rt_profiling.getProfilingInfo<CL_PROFILING_COMMAND_START>(&err);
                    CHECK_CL_ERR(err);
                    cl_ulong end = rt_profiling.getProfilingInfo<CL_PROFILING_COMMAND_END>(&err);
                    CHECK_CL_ERR(err);
                    duration = end - start;
                }
#else
                err = queue->enqueueNDRangeKernel(  refit_tree_kernel
                                                  , ::cl::NullRange
                                                  , ::cl::NDRange(total_objects * rt_local_work_size)
                                                  , ::cl::NDRange(rt_local_work_size)
                                                  );
                CHECK_CL_ERR(err);
#endif // USE_PROFILING
                return duration;
            }
            
        } // namespace cl
        
    } // namespace details
    
    template<typename V, size_t K, typename T>
    inline void refit_tree(  Tree<T,K> & tree
                             , mesh_array::T4Mesh const & mesh
                             , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & X
                             , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & Y
                             , mesh_array::VertexAttribute<T,mesh_array::T4Mesh> const & Z
                             , const dikucl & /* tag */
                             )
    {
        assert( 0 || !"refit_tree using DIKUCL is not yet implemented.");
        
        // TODO Batch refitting of trees, copy to device, call above routine
        // and copy back from device.
    }
    
} // namespace kdop

#endif // KDOP_CL_REFIT_TREE_H