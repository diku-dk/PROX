#ifndef __GLOBAL_EXACT_TESTS_CAPACITY
#error "__GLOBAL_EXACT_TESTS_CAPACITY must be defined"
#endif // __GLOBAL_EXACT_TESTS_CAPACITY

#ifndef __GLOBAL_WORK_QUEUE_CAPACITY
#error "__GLOBAL_WORK_QUEUE_CAPACITY must be defined"
#endif // __GLOBAL_WORK_QUEUE_CAPACITY

#ifndef __IDLE_THRESHOLD
#error "__IDLE_THRESHOLD must be defined"
#endif // __IDLE_THRESHOLD

#ifndef __LOCAL_WORK_QUEUE_CAPACITY
#error "__LOCAL_WORK_QUEUE_CAPACITY must be defined"
#endif // __LOCAL_WORK_QUEUE_CAPACITY

#ifndef __LOCAL_WORK_QUEUE_INIT_ITEMS
#error "__LOCAL_WORK_QUEUE_INIT_ITEMS must be defined"
#endif // __LOCAL_WORK_QUEUE_INIT_ITEMS

#ifndef __TRAVERSAL_THREADS
#error "__TRAVERSAL_THREADS must be defined"
#endif // __TRAVERSAL_THREADS

#include "kdop_cl_index_type.clh"
#include "kdop_cl_node.clh"
#include "kdop_cl_work_item.clh"

__kernel void do_tandem_traversal(
    __global const    Node     *restrict nodes,
    __global          WorkItem *restrict global_work_queues,
    __global          uint     *restrict global_work_queues_counts,
    __global volatile uint     *restrict idle_count,
    __global          WorkItem *restrict exact_tests,
    __global volatile uint     *restrict exact_tests_count)
{

    __local int local_work_queue_size;
    __local int global_work_queue_size;
    __local WorkItem local_work_queue[__LOCAL_WORK_QUEUE_CAPACITY];
    __local uint want_to_abort;

    const int group_id = get_group_id(0);
    const int local_id = get_local_id(0);

    if(local_id == 0) {
        global_work_queue_size = global_work_queues_counts[group_id];
        local_work_queue_size = min(__LOCAL_WORK_QUEUE_INIT_ITEMS,
                                    global_work_queue_size);
        global_work_queue_size -= local_work_queue_size;
    }
    
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

    if(local_work_queue_size == 0) {
        if(local_id == 0) {
            (void) atomic_inc(idle_count);
        }
        return;
    }

    {
        __global WorkItem *global_queue = &global_work_queues[group_id * __GLOBAL_WORK_QUEUE_CAPACITY + global_work_queue_size];
        int queue_offset = local_id;
        while(queue_offset < local_work_queue_size/* * 3*/) {
            local_work_queue[queue_offset] = global_queue[queue_offset];
            queue_offset += __TRAVERSAL_THREADS;
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    
    while(local_work_queue_size > 0) {
        WorkItem work_item;
        int n_active = min(__TRAVERSAL_THREADS, local_work_queue_size);

        bool can_work = local_id < local_work_queue_size;
        if(can_work) {
            work_item = local_work_queue[local_work_queue_size - n_active + local_id];
        }
        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

        if(local_id == 0) {
            local_work_queue_size -= n_active;
        }
        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
        
        if(can_work) {
            Node node_a = nodes[work_item.a];
            Node node_b = nodes[work_item.b];

            bool overlap = node_overlap(node_a, node_b);
            bool a_is_leaf = node_is_leaf(node_a);
            bool b_is_leaf = node_is_leaf(node_b);

            if(overlap) {
                if(a_is_leaf && b_is_leaf) {
                    uint exact_tests_pos = atomic_inc(exact_tests_count);
                    if(exact_tests_pos >= __GLOBAL_EXACT_TESTS_CAPACITY) {
                        // TODO handle
                    }
                    exact_tests[exact_tests_pos].a = node_a.start;
                    exact_tests[exact_tests_pos].b = node_b.start;
                    exact_tests[exact_tests_pos].tp = work_item.tp;
                } else if(b_is_leaf || (!a_is_leaf && node_size(node_a) > node_size(node_b))) {
                    int l = atomic_add(&local_work_queue_size, 2);

                    if(l < __LOCAL_WORK_QUEUE_CAPACITY) {
                        local_work_queue[l].a = node_a.start;
                        local_work_queue[l].b = work_item.b;
                        local_work_queue[l].tp = work_item.tp;
                        ++l;

                        if(l < __LOCAL_WORK_QUEUE_CAPACITY) {
                            local_work_queue[l].a = node_a.end;
                            local_work_queue[l].b = work_item.b;
                            local_work_queue[l].tp = work_item.tp;
                        } else {
                            int g = atomic_inc(&global_work_queue_size);
                            __global WorkItem *global_queue = &global_work_queues[group_id * __GLOBAL_WORK_QUEUE_CAPACITY + g];

                            global_queue[0].a = node_a.end;
                            global_queue[0].b = work_item.b;
                            global_queue[0].tp = work_item.tp;
                        }
                    } else {
                        int g = atomic_add(&global_work_queue_size, 2);
                        __global WorkItem *global_queue = &global_work_queues[group_id * __GLOBAL_WORK_QUEUE_CAPACITY + g];

                        global_queue[0].a = node_a.start;
                        global_queue[0].b = work_item.b;
                        global_queue[0].tp = work_item.tp;
                        global_queue[1].a = node_a.end;
                        global_queue[1].b = work_item.b;
                        global_queue[1].tp = work_item.tp;
                    }
                } else {
                    int l = atomic_add(&local_work_queue_size, 2);

                    if(l < __LOCAL_WORK_QUEUE_CAPACITY) {
                        local_work_queue[l].a = work_item.a;
                        local_work_queue[l].b = node_b.start;
                        local_work_queue[l].tp = work_item.tp;
                        ++l;

                        if(l < __LOCAL_WORK_QUEUE_CAPACITY) {
                            local_work_queue[l].a = work_item.a;
                            local_work_queue[l].b = node_b.end;
                            local_work_queue[l].tp = work_item.tp;
                        } else {
                            int g = atomic_inc(&global_work_queue_size);
                            __global WorkItem *global_queue = &global_work_queues[group_id * __GLOBAL_WORK_QUEUE_CAPACITY + g];

                            global_queue[0].a = work_item.a;
                            global_queue[0].b = node_b.end;
                            global_queue[0].tp = work_item.tp;
                        }
                    } else {
                        int g = atomic_add(&global_work_queue_size, 2);
                        __global WorkItem *global_queue = &global_work_queues[group_id * __GLOBAL_WORK_QUEUE_CAPACITY + g];

                        global_queue[0].a = work_item.a;
                        global_queue[0].b = node_b.start;
                        global_queue[0].tp = work_item.tp;
                        global_queue[1].a = work_item.a;
                        global_queue[1].b = node_b.end;
                        global_queue[1].tp = work_item.tp;
                    }
                }
            }
        }
        
        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

        if(local_work_queue_size >= __LOCAL_WORK_QUEUE_CAPACITY - __TRAVERSAL_THREADS ||
           global_work_queue_size >= __GLOBAL_WORK_QUEUE_CAPACITY - __TRAVERSAL_THREADS * 2 - __LOCAL_WORK_QUEUE_CAPACITY ||
           local_work_queue_size == 0)
        {
            if(local_id == 0) {
                (void) atomic_inc(idle_count);
            }
            break;
        }

        if(local_id == 0) {
            want_to_abort = *idle_count;
        }

        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

        if(want_to_abort >= __IDLE_THRESHOLD) {
            if(local_id == 0) {
                (void) atomic_inc(idle_count);
            }
            break;
        }
    }
    
    if(local_id == 0) {
        local_work_queue_size = min(local_work_queue_size, __LOCAL_WORK_QUEUE_CAPACITY);
        global_work_queues_counts[group_id] = local_work_queue_size + global_work_queue_size;
    }
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

    {
        int queue_offset = local_id;
        __global WorkItem *global_queue = &global_work_queues[group_id * __GLOBAL_WORK_QUEUE_CAPACITY + global_work_queue_size];
        while(queue_offset < local_work_queue_size/* * 2*/) {
            global_queue[queue_offset] = local_work_queue[queue_offset];
            queue_offset += __TRAVERSAL_THREADS;
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
}
// prevent file parsing warning on AMD APP SDK v2.9 by appending this line