.. _interface_TKE:

TKE
===

Inside the library, the Fortran derived data types are rebuilt as lightweight structures. This is done to improve the code readability and to simplify the interface to the library backend.

For example, internally the library defines a `t_patch` structures as follows::

   struct t_patch {
       double *depth_CellInterface;
       double *prism_center_dist_c;
       double *inv_prism_center_dist_c;
       double *prism_thick_c;
       int *dolic_c;
       int *dolic_e;
       double *zlev_i;
       double *wet_c;
       int *edges_cell_idx;
       int *edges_cell_blk;
   };

which is a lightweight version of the ICON version.

This structures are filled inside the `TKE` class during the first time step. Then, these structures are provided to the backend.

The backend is internally using a memory view on the allocated memory. The interface allows to use different memory views implementations for different backends (CPU, CUDA or HIP) and to easily change to a different memory view implementation from an existing one. For example, the CUDA backend uses the `mdspan` from the CUDA standard library which can generate a 1D, 2D or 3D view based on a provided memory allocation and it allows to use the allocated contiguous one dimensional memory as Fortran arrays. The memory view objects are created during the first time step based on the pointers provided by the model and they are organized in structures of memory views. These structures are then used in the computations. 

In order to achieve enough flexibility in the interface, the structures of memory views are templated. For example a structure of memory views mirroring the `t_patch` struct::

   template <template <class ...> class memview,
             template <class, size_t> class dext>
   struct t_patch_view {
       memview<double, dext<int, 3>> depth_CellInterface;
       memview<double, dext<int, 3>> prism_center_dist_c;
       memview<double, dext<int, 3>> inv_prism_center_dist_c;
       memview<double, dext<int, 3>> prism_thick_c;
       memview<int, dext<int, 2>> dolic_c;
       memview<int, dext<int, 2>> dolic_e;
       memview<double, dext<int, 1>> zlev_i;
       memview<double, dext<int, 3>> wet_c;
       memview<int, dext<int, 3>> edges_cell_idx;
       memview<int, dext<int, 3>> edges_cell_blk;
   };

The assumption is that a memory view class is defined by a type and an extent class. It is then possible to write an implementation of those two classes and use it in the TKE library interface.

The main benefit of the memory views approach is that the actual code looks very similar to the original ICON code written in Fortran. For example, an iteration over `prism_thick_c`::

    for (int jb = cells_start_block; jb <= cells_end_block; jb++) {
        int start_index, end_index;
        get_index_range(cells_block_size, cells_start_block, cells_end_block,
                        cells_start_index, cells_end_index, jb, &start_index, &end_index);
        for (int level = 0; level < nlevels; level++) {
            for (int jc = start_index; jc <= end_index; jc++) {
                p_patch.prism_thick_c(jb, level, jc) = ...
            }
        }
    }

In order to be able to use only memory views inside the vertical mixing scheme, the internal fields are allocated during the initilization step and then a mdspan object is created based on the allocated memory pointer.

Each backend defines a memory view policy class which defines the methods to create a memory view object given a pointer and to allocate memory and associate a memory view object to it (policy based design).
