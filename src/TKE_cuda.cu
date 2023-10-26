/* Copyright (C) 2023  Enrico Degregori, Wilton Jaciel Loch
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "src/TKE_cuda.hpp"


#include <cuda.h>
#include <cuda/std/mdspan>
#include <iostream>

#include "src/utils.hpp"
#include "src/cuda_check.hpp"

constexpr auto dyn = cuda::std::dynamic_extent;
using ext1d_t = cuda::std::extents<size_t, dyn>;
using ext2d_t = cuda::std::extents<size_t, dyn, dyn>;
using ext3d_t = cuda::std::extents<size_t, dyn, dyn, dyn>;
using mdspan_1d_double = cuda::std::mdspan<double, ext1d_t>;
using mdspan_2d_double = cuda::std::mdspan<double, ext2d_t>;
using mdspan_3d_double = cuda::std::mdspan<double, ext3d_t>;
using mdspan_2d_int = cuda::std::mdspan<int, ext2d_t>;

// this needs to be defined static or TKE_cuda.hpp can not
// be included in C++ files (need to think about other solutions)
static mdspan_1d_double view_cuda_malloc(double *field, size_t dim1);
static mdspan_2d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2);

// TKE internal memory views
static mdspan_2d_double rho_up_view;
static mdspan_2d_double rho_down_view;
static mdspan_1d_double forc_tke_surf_2D_view;
static mdspan_1d_double forc_rho_surf_2D_view;
static mdspan_1d_double bottom_fric_2D_view;
static mdspan_1d_double s_c_view;
static mdspan_2d_double dzw_stretched_view;
static mdspan_2d_double dzt_stretched_view;
static mdspan_2d_double tke_old_view;
static mdspan_2d_double tke_Av_view;
static mdspan_2d_double tke_kv_view;
static mdspan_2d_double tke_iw_alpha_c_view;
static mdspan_2d_double tke_iwe_view;
static mdspan_2d_double tke_iwe_forcing_view;
static mdspan_2d_double pressure_view;
static mdspan_2d_double Nsqr_view;
static mdspan_2d_double Ssqr_view;

// TKE interface memory views
static mdspan_3d_double tke_view;
static mdspan_2d_int dolic_c_view;

struct t_cvmix_view {
    mdspan_3d_double tke;
    mdspan_3d_double tke_plc;
    mdspan_2d_double hlc;
    mdspan_3d_double wlc;
    mdspan_2d_double u_stokes;
    mdspan_3d_double a_veloc_v;
    mdspan_3d_double a_temp_v;
    mdspan_3d_double a_salt_v;
    mdspan_3d_double iwe_Tdis;
    mdspan_3d_double cvmix_dummy_1;
    mdspan_3d_double cvmix_dummy_2;
    mdspan_3d_double cvmix_dummy_3;
    mdspan_3d_double tke_Tbpr;
    mdspan_3d_double tke_Tspr;
    mdspan_3d_double tke_Tdif;
    mdspan_3d_double tke_Tdis;
    mdspan_3d_double tke_Twin;
    mdspan_3d_double tke_Tiwf;
    mdspan_3d_double tke_Tbck;
    mdspan_3d_double tke_Ttot;
    mdspan_3d_double tke_Lmix;
    mdspan_3d_double tke_Pr;
};

struct t_cvmix_view *p_cvmix_view;

static void fill_struct_view(struct t_cvmix_view *p_cvmix_view_d, struct t_cvmix *p_cvmix,
                             int nblocks, int nlevs, int nproma);

// TKE CUDA kernels functions
__global__ void calc_impl_kernel(int blockNo, int start_index, int end_index,
                                 mdspan_2d_int dolic_c, struct t_cvmix_view *p_cvmix,
                                 mdspan_2d_double tke_old);

TKE_cuda::TKE_cuda(int nproma, int nlevs, int nblocks, int vert_mix_type, int vmix_idemix_tke,
                   int vert_cor_type, double dtime, double OceanReferenceDensity, double grav,
                   int l_lc, double clc, double ReferencePressureIndbars, double pi)
    : TKE_backend(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke,
                  vert_cor_type, dtime, OceanReferenceDensity, grav,
                  l_lc, clc, ReferencePressureIndbars, pi) {
    // Initialize internal arrays
    std::cout << "Initializing TKE cuda... " << std::endl;
    rho_up_view = view_cuda_malloc(m_rho_up, static_cast<size_t>(nlevs), static_cast<size_t>(nproma));
    rho_up_view = view_cuda_malloc(m_rho_down, static_cast<size_t>(nlevs), static_cast<size_t>(nproma));
    forc_tke_surf_2D_view = view_cuda_malloc(m_forc_tke_surf_2D, static_cast<size_t>(nproma));
    forc_rho_surf_2D_view = view_cuda_malloc(m_forc_rho_surf_2D, static_cast<size_t>(nproma));
    bottom_fric_2D_view = view_cuda_malloc(m_bottom_fric_2D, static_cast<size_t>(nproma));
    s_c_view = view_cuda_malloc(m_s_c, static_cast<size_t>(nproma));
    dzw_stretched_view = view_cuda_malloc(m_dzw_stretched, static_cast<size_t>(nlevs), static_cast<size_t>(nproma));
    dzt_stretched_view = view_cuda_malloc(m_dzt_stretched, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    tke_old_view = view_cuda_malloc(m_tke_old, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    tke_Av_view = view_cuda_malloc(m_tke_Av, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    tke_kv_view = view_cuda_malloc(m_tke_kv, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    tke_iw_alpha_c_view = view_cuda_malloc(m_tke_iw_alpha_c, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    tke_iwe_view = view_cuda_malloc(m_tke_iwe, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    tke_iwe_forcing_view = view_cuda_malloc(m_tke_iwe_forcing, static_cast<size_t>(nlevs+1),
                                            static_cast<size_t>(nproma));
    pressure_view = view_cuda_malloc(m_pressure, static_cast<size_t>(nlevs), static_cast<size_t>(nproma));
    Nsqr_view = view_cuda_malloc(m_Nsqr, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));
    Ssqr_view = view_cuda_malloc(m_Ssqr, static_cast<size_t>(nlevs+1), static_cast<size_t>(nproma));

    is_view_init = false;
}

TKE_cuda::~TKE_cuda() {
    // Free internal arrays memory
    std::cout << "Finalizing TKE cuda... " << std::endl;
    check(cudaFree(m_rho_up) );
    check(cudaFree(m_rho_down) );
    check(cudaFree(m_forc_tke_surf_2D) );
    check(cudaFree(m_forc_rho_surf_2D) );
    check(cudaFree(m_bottom_fric_2D) );
    check(cudaFree(m_s_c) );
    check(cudaFree(m_dzw_stretched) );
    check(cudaFree(m_dzt_stretched) );
    check(cudaFree(m_tke_old) );
    check(cudaFree(m_tke_Av) );
    check(cudaFree(m_tke_kv) );
    check(cudaFree(m_tke_iw_alpha_c) );
    check(cudaFree(m_tke_iwe) );
    check(cudaFree(m_tke_iwe_forcing) );
    check(cudaFree(m_pressure) );
    check(cudaFree(m_Nsqr) );
    check(cudaFree(m_Ssqr) );
}

void TKE_cuda::calc_impl(struct t_patch p_patch, struct t_cvmix p_cvmix,
                         struct t_ocean_state ocean_state, struct t_atmo_fluxes atmos_fluxes,
                         struct t_atmos_for_ocean p_as, struct t_sea_ice p_sea_ice,
                         int edges_block_size, int edges_start_block, int edges_end_block,
                         int edges_start_index, int edges_end_index, int cells_block_size,
                         int cells_start_block, int cells_end_block, int cells_start_index,
                         int cells_end_index) {
    if (!is_view_init) {
        fill_struct_view(p_cvmix_view, &p_cvmix, m_nblocks, m_nlevs, m_nproma);
        dolic_c_view = mdspan_2d_int{ p_patch.dolic_c, ext2d_t{m_nblocks, m_nproma} };
        is_view_init = true;
    }

    for (int jb = cells_start_block; jb <= cells_end_block; jb++) {
        int start_index, end_index;
        get_index_range(cells_block_size, cells_start_block, cells_end_block,
                        cells_start_index, cells_end_index, jb, &start_index, &end_index);
        int threadsPerBlockI = 512;
        int blocksPerGridI = (end_index - start_index) / threadsPerBlockI + 1;
        dim3 blocksPerGrid(blocksPerGridI, 1, 1);
        dim3 threadsPerBlock(threadsPerBlockI, 1, 1);
        calc_impl_kernel<<<blocksPerGrid, threadsPerBlock>>>(jb, start_index, end_index,
                                                            dolic_c_view, p_cvmix_view,
                                                            tke_old_view);
    }
                         }

static mdspan_1d_double view_cuda_malloc(double *field, size_t dim1) {
    check( cudaMalloc(&field, dim1*sizeof(double)) );
    mdspan_1d_double memview{ field, ext1d_t{dim1} };
    return memview;
}

static mdspan_2d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2) {
    check( cudaMalloc(&field, dim1*dim2*sizeof(double)) );
    mdspan_2d_double memview{ field, ext2d_t{dim1, dim2} };
    return memview;
}

static void fill_struct_view(struct t_cvmix_view *p_cvmix_view_d, struct t_cvmix *p_cvmix,
                             int nblocks, int nlevs, int nproma) {
    struct t_cvmix_view p_cvmix_view_l;
    // create mdspan CPU object containing GPU pointers provided by the frontend
    p_cvmix_view_l.tke = mdspan_3d_double{ p_cvmix->tke, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_plc = mdspan_3d_double{ p_cvmix->tke_plc, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.hlc = mdspan_2d_double{ p_cvmix->hlc, ext2d_t{nblocks, nproma} };
    p_cvmix_view_l.wlc = mdspan_3d_double{ p_cvmix->wlc, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.u_stokes = mdspan_2d_double{ p_cvmix->u_stokes, ext2d_t{nblocks, nproma} };
    p_cvmix_view_l.a_veloc_v = mdspan_3d_double{ p_cvmix->a_veloc_v, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.a_temp_v = mdspan_3d_double{ p_cvmix->a_temp_v, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.a_salt_v = mdspan_3d_double{ p_cvmix->a_salt_v, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.iwe_Tdis = mdspan_3d_double{ p_cvmix->iwe_Tdis, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.cvmix_dummy_1 = mdspan_3d_double{ p_cvmix->cvmix_dummy_1, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.cvmix_dummy_2 = mdspan_3d_double{ p_cvmix->cvmix_dummy_2, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.cvmix_dummy_3 = mdspan_3d_double{ p_cvmix->cvmix_dummy_3, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Tbpr = mdspan_3d_double{ p_cvmix->tke_Tbpr, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Tspr = mdspan_3d_double{ p_cvmix->tke_Tspr, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Tdif = mdspan_3d_double{ p_cvmix->tke_Tdif, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Tdis = mdspan_3d_double{ p_cvmix->tke_Tdis, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Twin = mdspan_3d_double{ p_cvmix->tke_Twin, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Tiwf = mdspan_3d_double{ p_cvmix->tke_Tiwf, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Tbck = mdspan_3d_double{ p_cvmix->tke_Tbck, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Ttot = mdspan_3d_double{ p_cvmix->tke_Ttot, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Lmix = mdspan_3d_double{ p_cvmix->tke_Lmix, ext3d_t{nblocks, nlevs, nproma} };
    p_cvmix_view_l.tke_Pr = mdspan_3d_double{ p_cvmix->tke_Pr, ext3d_t{nblocks, nlevs, nproma} };
    // allocate t_cvmix_view structure on the GPU
    check( cudaMalloc(&p_cvmix_view_d, sizeof(t_cvmix_view)) );
    // copy CPU t_cvmix_view to GPU
    check( cudaMemcpy(p_cvmix_view_d, &p_cvmix_view_l, sizeof(t_cvmix_view), cudaMemcpyHostToDevice) );
}

__global__ void calc_impl_kernel(int blockNo, int start_index, int end_index,
                                 mdspan_2d_int dolic_c, struct t_cvmix_view *p_cvmix,
                                 mdspan_2d_double tke_old) {
    int jc = blockIdx.x * blockDim.x + threadIdx.x + start_index;
    if (jc <= end_index) {
        int levels = dolic_c(blockNo, jc);
        for (int level = 0; level < levels; level++) {
            tke_old(level, jc) = p_cvmix->tke(blockNo, level, jc);
            p_cvmix->tke(blockNo, level, jc) = p_cvmix->tke(blockNo, level, jc) + 1.0;
        }
    }
}
