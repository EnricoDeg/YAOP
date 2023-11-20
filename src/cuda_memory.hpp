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

#ifndef SRC_CUDA_MEMORY_HPP_
#define SRC_CUDA_MEMORY_HPP_

#include <cuda.h>
#include <cuda/std/mdspan>
#include "src/cuda_check.hpp"
#include "src/data_struct.hpp"

constexpr auto dyn = cuda::std::dynamic_extent;
using ext1d_t = cuda::std::extents<size_t, dyn>;
using ext2d_t = cuda::std::extents<size_t, dyn, dyn>;
using ext3d_t = cuda::std::extents<size_t, dyn, dyn, dyn>;
using mdspan_1d_double = cuda::std::mdspan<double, ext1d_t>;
using mdspan_2d_double = cuda::std::mdspan<double, ext2d_t>;
using mdspan_3d_double = cuda::std::mdspan<double, ext3d_t>;
using mdspan_2d_int = cuda::std::mdspan<int, ext2d_t>;
using mdspan_3d_int = cuda::std::mdspan<int, ext3d_t>;

// TKE constants
struct t_constant {
    int vert_mix_type;
    int vmix_idemix_tke;
    int vert_cor_type;
    double dtime;
    double OceanReferenceDensity;
    double grav;
    int l_lc;
    double clc;
    double ReferencePressureIndbars;
    double pi;
    int nlevs;
};

struct t_constant_tke {
    double c_k;
    double c_eps;
    double cd;
    double alpha_tke;
    double clc;
    double mxl_min;
    double KappaM_min;
    double KappaH_min;
    double KappaM_max;
    double tke_surf_min;
    double tke_min;
    int tke_mxl_choice;
    int handle_old_vals;
    bool only_tke;
    bool use_Kappa_min;
    bool use_ubound_dirichlet;
    bool use_lbound_dirichlet;
};

// TKE internal memory views
struct t_tke_internal_view {
    mdspan_2d_double rho_up;
    mdspan_2d_double rho_down;
    mdspan_1d_double forc_tke_surf_2D;
    mdspan_1d_double forc_rho_surf_2D;
    mdspan_1d_double bottom_fric_2D;
    mdspan_1d_double s_c;
    mdspan_2d_double dzw_stretched;
    mdspan_2d_double dzt_stretched;
    mdspan_2d_double tke_old;
    mdspan_3d_double tke_Av;
    mdspan_2d_double tke_kv;
    mdspan_2d_double tke_iwe;
    mdspan_2d_double tke_iwe_forcing;
    mdspan_2d_double pressure;
    mdspan_2d_double Nsqr;
    mdspan_2d_double Ssqr;
    mdspan_2d_double a_dif;
    mdspan_2d_double b_dif;
    mdspan_2d_double c_dif;
    mdspan_2d_double a_tri;
    mdspan_2d_double b_tri;
    mdspan_2d_double c_tri;
    mdspan_2d_double d_tri;
    mdspan_2d_double mxl;
    mdspan_2d_double sqrttke;
    mdspan_2d_double KappaM_out;
    mdspan_2d_double Rinum;
    mdspan_2d_double prandtl;
    mdspan_2d_double KappaH_out;
    mdspan_2d_double forc;
    mdspan_2d_double K_diss_v;
    mdspan_2d_double P_diss_v;
    mdspan_2d_double ke;
    mdspan_2d_double cp;
    mdspan_2d_double dp;
    mdspan_2d_double tke_upd;
    mdspan_2d_double tke_unrest;
};

// TKE interface memory views
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

struct t_patch_view {
    mdspan_3d_double depth_CellInterface;
    mdspan_3d_double prism_center_dist_c;
    mdspan_3d_double inv_prism_center_dist_c;
    mdspan_3d_double prism_thick_c;
    mdspan_2d_int dolic_c;
    mdspan_2d_int dolic_e;
    mdspan_1d_double zlev_i;
    mdspan_3d_double wet_c;
    mdspan_3d_int edges_cell_idx;
    mdspan_3d_int edges_cell_blk;
};

struct t_ocean_state_view {
    mdspan_3d_double temp;
    mdspan_3d_double salt;
    mdspan_2d_double stretch_c;
    mdspan_2d_double eta_c;
    mdspan_3d_double p_vn_x1;
    mdspan_3d_double p_vn_x2;
    mdspan_3d_double p_vn_x3;
};

struct t_atmo_fluxes_view {
    mdspan_2d_double stress_xw;
    mdspan_2d_double stress_yw;
};

struct t_atmos_for_ocean_view {
    mdspan_2d_double fu10;
};

struct t_sea_ice_view {
    mdspan_2d_double concsum;
};

void fill_struct_view(struct t_cvmix_view *p_cvmix_view, struct t_cvmix *p_cvmix,
                             int nblocks, int nlevs, int nproma);

void fill_struct_view(struct t_patch_view *p_patch_view, struct t_patch *p_patch,
                             int nblocks, int nlevs, int nproma);


void fill_struct_view(struct t_ocean_state_view *ocean_state_view, struct t_ocean_state *ocean_state,
                             int nblocks, int nlevs, int nproma);

void fill_struct_view(struct t_atmo_fluxes_view *atmos_fluxes_view, struct t_atmo_fluxes *atmos_fluxes,
                             int nblocks, int nlevs, int nproma);

void fill_struct_view(struct t_atmos_for_ocean_view *p_as_view, struct t_atmos_for_ocean *p_as,
                             int nblocks, int nlevs, int nproma);

void fill_struct_view(struct t_sea_ice_view *p_sea_ice_view, struct t_sea_ice *p_sea_ice,
                             int nblocks, int nlevs, int nproma);

mdspan_1d_double view_cuda_malloc(double *field, size_t dim1);

mdspan_2d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2);

mdspan_3d_double view_cuda_malloc(double *field, size_t dim1, size_t dim2, size_t dim3);

#endif  // SRC_CUDA_MEMORY_HPP_
