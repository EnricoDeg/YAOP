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

#ifndef SRC_BACKENDS_CPU_CPU_KERNELS_HPP_
#define SRC_BACKENDS_CPU_CPU_KERNELS_HPP_

#include <algorithm>
#include <cmath>
#include "src/backends/CPU/cpu_memory.hpp"
#include "src/shared/interface/memview_struct.hpp"
#include "src/shared/interface/data_struct.hpp"
#include "src/backends/kernels.hpp"
#include "src/shared/constants/constants_thermodyn.hpp"

using std::max;
using std::min;

template <class T>
void calc_impl_edges(int blockNo, int start_index, int end_index,
                     t_patch_view<T, cpu_memview::mdspan, cpu_memview::dextents> p_patch,
                     t_cvmix_view<T, cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
                     t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
                     t_constant p_constant) {
    // compute max level on block (maxval fortran function)
    int max_levels = 0;
    for (int je = start_index; je <= end_index; je++)
        if (p_patch.dolic_e(blockNo, je) > max_levels)
            max_levels = p_patch.dolic_e(blockNo, je);

    for (int level = 1; level < p_constant.nlevs+1; level++) {
        for (int je = start_index; je <= end_index; je++) {
            int cell_1_idx = p_patch.edges_cell_idx(0, blockNo, je);
            int cell_1_block = p_patch.edges_cell_blk(0, blockNo, je);
            int cell_2_idx = p_patch.edges_cell_idx(1, blockNo, je);
            int cell_2_block = p_patch.edges_cell_blk(1, blockNo, je);

            if (level < p_patch.dolic_e(blockNo, je))
                p_cvmix.a_veloc_v(blockNo, level, je) = 0.5 * (p_internal.tke_Av(cell_1_block, level, cell_1_idx) +
                                                        p_internal.tke_Av(cell_2_block, level, cell_2_idx));
            else
                p_cvmix.a_veloc_v(blockNo, level, je) = 0.0;
        }
    }
}

template <class T>
inline
void calc_mxl_2(int blockNo, int start_index, int end_index, int max_levels, T mxl_min,
                mdspan_2d<int> dolic_c, mdspan_3d<T> tke_Lmix, mdspan_2d<T> dzw_stretched) {
    for (int jc = start_index; jc <= end_index; jc++) {
        if (dolic_c(blockNo, jc) > 0) {
            tke_Lmix(blockNo, 0, jc) = 0.0;
            tke_Lmix(blockNo, dolic_c(blockNo, jc), jc) = 0.0;
        }
    }

    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                tke_Lmix(blockNo, level, jc) = min(tke_Lmix(blockNo, level, jc),
                         tke_Lmix(blockNo, level-1, jc) + dzw_stretched(level-1, jc));

    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc) > 0) {
            int dolic = dolic_c(blockNo, jc);
            tke_Lmix(blockNo, dolic-1, jc) = min(tke_Lmix(blockNo, dolic-1, jc),
                                             mxl_min + dzw_stretched(dolic-1, jc));
        }

    for (int level = max_levels-2; level > 0; level--)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc) - 1)
                tke_Lmix(blockNo, level, jc) = min(tke_Lmix(blockNo, level, jc),
                         tke_Lmix(blockNo, level+1, jc) +  dzw_stretched(level, jc));

    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc) + 1)
                tke_Lmix(blockNo, level, jc) = max(tke_Lmix(blockNo, level, jc), mxl_min);
}

template <class T>
inline
void calc_diffusivity(int blockNo, int start_index, int end_index, int max_levels,
                      t_constant_tke *p_constant_tke,
                      mdspan_2d<int> dolic_c, mdspan_3d<T> tke_Lmix, mdspan_2d<T> sqrttke,
                      mdspan_2d<T> Nsqr, mdspan_2d<T> Ssqr,
                      mdspan_3d<T> tke_Av, mdspan_2d<T> tke_kv, mdspan_3d<T> tke_Pr) {
    for (int level = 0; level < max_levels+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < dolic_c(blockNo, jc) + 1) {
                tke_Av(blockNo, level, jc) = min(p_constant_tke->KappaM_max,
                                                 p_constant_tke->c_k * tke_Lmix(blockNo, level, jc) *
                                                 sqrttke(level, jc));
                tke_Pr(blockNo, level, jc) = Nsqr(level, jc) / max(Ssqr(level, jc), 1.0e-12);
                if (!p_constant_tke->only_tke)
                    tke_Pr(blockNo, level, jc) = min(tke_Pr(blockNo, level, jc),
                                                     tke_Av(blockNo, level, jc) * Nsqr(level, jc) / 1.0e-12);
                tke_Pr(blockNo, level, jc) = max(1.0, min(10.0, 6.6 * tke_Pr(blockNo, level, jc)));
                tke_kv(level, jc) = tke_Av(blockNo, level, jc) / tke_Pr(blockNo, level, jc);
                if (p_constant_tke->use_Kappa_min) {
                    tke_Av(blockNo, level, jc) = max(p_constant_tke->KappaM_min, tke_Av(blockNo, level, jc));
                    tke_kv(level, jc) = max(p_constant_tke->KappaH_min, tke_kv(level, jc));
                }
            }
        }
    }
}

template <class T>
inline
void calc_forcing(int blockNo, int start_index, int end_index, int max_levels, bool l_lc, bool only_tke,
                  mdspan_2d<int> dolic_c, mdspan_2d<T> Ssqr, mdspan_2d<T> Nsqr, mdspan_3d<T> tke_Av,
                  mdspan_2d<T> tke_kv, mdspan_3d<T> tke_Tspr, mdspan_3d<T> tke_Tbpr,
                  mdspan_3d<T> tke_plc, mdspan_3d<T> tke_Tiwf, mdspan_2d<T> forc) {
    for (int level = 0; level < max_levels+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < dolic_c(blockNo, jc) + 1) {
                // forcing by shear and buoycancy production
                tke_Tspr(blockNo, level, jc) = Ssqr(level, jc) * tke_Av(blockNo, level, jc);
                tke_Tbpr(blockNo, level, jc) = Nsqr(level, jc) * tke_kv(level, jc);
                if (level == 0) tke_Tbpr(blockNo, 0, jc) = 0.0;

                forc(level, jc) = tke_Tspr(blockNo, level, jc) - tke_Tbpr(blockNo, level, jc);
                // additional langmuir turbulence term
                if (l_lc)
                    forc(level, jc) += tke_plc(blockNo, level, jc);
                // forcing by internal wave dissipation
                if (!only_tke)
                    forc(level, jc) += tke_Tiwf(blockNo, level, jc);
            }
        }
    }
}

template <class T>
inline
void build_diffusion_dissipation_tridiag(int blockNo, int start_index, int end_index, int max_levels,
                                         mdspan_2d<int> dolic_c, T alpha_tke,
                                         mdspan_3d<T> tke_Av, mdspan_2d<T> dzt_stretched,
                                         mdspan_2d<T> dzw_stretched,
                                         mdspan_2d<T> ke, mdspan_2d<T> a_dif, mdspan_2d<T> b_dif,
                                         mdspan_2d<T> c_dif) {
    // c is lower diagonal of matrix
    for (int level = 0; level < max_levels; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < dolic_c(blockNo, jc)) {
                int kp1 = min(level+1, dolic_c(blockNo, jc)-1);
                int kk = max(level, 1);
                ke(level, jc) = 0.5 * alpha_tke * (tke_Av(blockNo, kp1, jc) + tke_Av(blockNo, kk, jc));
                c_dif(level, jc) = ke(level, jc) / (dzt_stretched(level, jc) * dzw_stretched(level, jc));
            }
        }
    }

    // not part of the diffusion matrix, thus value is arbitrary (set to zero)
    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc) >= 0)
            c_dif(dolic_c(blockNo, jc), jc) = 0.0;

    // b is main diagonal of matrix
    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                b_dif(level, jc) = ke(level-1, jc) / (dzt_stretched(level, jc) * dzw_stretched(level-1, jc)) +
                                   ke(level, jc) / (dzt_stretched(level, jc) * dzw_stretched(level, jc));

    // a is upper diagonal of matrix
    for (int level = 1; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc) + 1)
                a_dif(level, jc) = ke(level-1, jc) / (dzt_stretched(level, jc) * dzw_stretched(level-1, jc));

    // not part of the diffusion matrix, thus value is arbitrary (set to zero)
    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc) >= 0)
            a_dif(0, jc) = 0.0;
}

template <class T>
inline
void build_tridiag(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                   T dtime, T c_eps, int nlevs,
                   mdspan_2d<T> a_dif, mdspan_2d<T> b_dif, mdspan_2d<T> c_dif,
                   mdspan_2d<T> sqrttke, mdspan_3d<T> tke_Lmix, mdspan_2d<T> tke_upd,
                   mdspan_2d<T> forc,
                   mdspan_2d<T> a_tri, mdspan_2d<T> b_tri, mdspan_2d<T> c_tri,
                   mdspan_2d<T> d_tri) {
    for (int level = 0; level < nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            a_tri(level, jc) = - dtime * a_dif(level, jc);
            b_tri(level, jc) = 1.0 + dtime * b_dif(level, jc);
            c_tri(level, jc) = - dtime * c_dif(level, jc);
        }
    }

    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                b_tri(level, jc) = b_tri(level, jc) + dtime * c_eps * sqrttke(level, jc) /
                                   tke_Lmix(blockNo, level, jc);

    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc) + 1)
                d_tri(level, jc) = tke_upd(level, jc) + dtime * forc(level, jc);
}

template <class T>
inline
void solve_tridiag(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                   mdspan_2d<T> a, mdspan_2d<T> b, mdspan_2d<T> c, mdspan_2d<T> d,
                   mdspan_3d<T> x, mdspan_2d<T> cp, mdspan_2d<T> dp) {
    // initialize a-prime (cp) and d-prime
    for (int jc = start_index; jc <= end_index; jc++) {
        if (dolic_c(blockNo, jc)+1 > 0) {
            int dolic = dolic_c(blockNo, jc);
            cp(dolic, jc) = c(dolic, jc) / b(dolic, jc);
            dp(dolic, jc) = d(dolic, jc) / b(dolic, jc);
        }
    }

    // solve for vectors a-prime and d-prime
    for (int level = max_levels-1; level >= 0; level--) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level <= dolic_c(blockNo, jc)) {
                T fxa = 1.0 / (b(level, jc) - cp(level+1, jc) * c(level, jc));
            }
        }
    }

    // initialize x
    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc)+1 > 0)
            x(blockNo, 0, jc) = dp(0, jc);

    // solve for x from the vectors a-prime and d-prime
    for (int level = 1; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level <= dolic_c(blockNo, jc))
                x(blockNo, level, jc) = dp(level, jc) - cp(level, jc) * x(blockNo, level-1, jc);
}

template <class T>
inline
void tke_vertical_diffusion(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                            T diff_surf_forc, T diff_bott_forc,
                            mdspan_2d<T> a_dif, mdspan_2d<T> b_dif, mdspan_2d<T> c_dif,
                            mdspan_3d<T> tke, mdspan_3d<T> tke_Tdif) {
    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                 tke_Tdif(blockNo, level, jc) = a_dif(level, jc) * tke(blockNo, level-1, jc) -
                                                b_dif(level, jc) * tke(blockNo, level, jc) +
                                                c_dif(level, jc) * tke(blockNo, level+1, jc);

    for (int jc = start_index; jc <= end_index; jc++) {
        if (dolic_c(blockNo, jc) > 0) {
            int dolic = dolic_c(blockNo, jc);
            tke_Tdif(blockNo, 0, jc) = - b_dif(0, jc) * tke(blockNo, 0, jc) +
                                         c_dif(0, jc) * tke(blockNo, 1, jc);
            tke_Tdif(blockNo, 1, jc) += diff_surf_forc;
            tke_Tdif(blockNo, dolic-1, jc) += diff_bott_forc;
            tke_Tdif(blockNo, dolic, jc) = a_dif(dolic, jc) * tke(blockNo, dolic-1, jc) -
                                           b_dif(dolic, jc) * tke(blockNo, dolic, jc);
        }
    }
}

template <class T>
inline
void tke_vertical_diffusion_ub_dirichlet(int blockNo, int start_index, int end_index, mdspan_2d<int> dolic_c,
                                         T tke_surf, mdspan_2d<T> ke, mdspan_2d<T> dzw_stretched,
                                         mdspan_2d<T> dzt_stretched, mdspan_3d<T> tke,
                                         mdspan_3d<T> tke_Tdif) {
    for (int jc = start_index; jc <= end_index; jc++)
        if (dolic_c(blockNo, jc) > 0)
            tke_Tdif(blockNo, 0, jc) = - ke(0, jc) / dzw_stretched(0, jc) / dzt_stretched(0, jc) *
                                       (tke_surf - tke(blockNo, 1, jc));
}

template <class T>
inline
void tke_vertical_diffusion_lb_dirichlet(int blockNo, int start_index, int end_index, mdspan_2d<int> dolic_c,
                                         T tke_bott, mdspan_2d<T> ke, mdspan_2d<T> dzw_stretched,
                                         mdspan_2d<T> dzt_stretched, mdspan_3d<T> tke,
                                         mdspan_3d<T> tke_Tdif) {
    for (int jc = start_index; jc <= end_index; jc++) {
        if (dolic_c(blockNo, jc) > 0) {
            int dolic = dolic_c(blockNo, jc);
            tke_Tdif(blockNo, dolic, jc) = ke(dolic-1, jc) / dzw_stretched(dolic-1, jc) /
                                           dzt_stretched(dolic, jc) * (tke(blockNo, dolic-1, jc) - tke_bott);
        }
    }
}

template <class T>
inline
void tke_vertical_dissipation(int blockNo, int start_index, int end_index, int max_levels, mdspan_2d<int> dolic_c,
                              int nlevs, T c_eps, mdspan_3d<T> tke_Lmix, mdspan_2d<T> sqrttke,
                              mdspan_3d<T> tke, mdspan_3d<T> tke_Tdis) {
    for (int level = 0; level < nlevs+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            tke_Tdis(blockNo, level, jc) = 0.0;

    for (int level = 1; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            if (level < dolic_c(blockNo, jc))
                tke_Tdis(blockNo, level, jc) = - c_eps / tke_Lmix(blockNo, level, jc) *
                                                 sqrttke(level, jc) * tke(blockNo, level, jc);
}

template <class T>
void integrate(int blockNo, int start_index, int end_index,
               t_patch_view<T, cpu_memview::mdspan, cpu_memview::dextents> p_patch,
               t_cvmix_view<T, cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
               t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
               t_constant p_constant,
               t_constant_tke p_constant_tke) {
    T tke_surf, diff_surf_forc, tke_bott, diff_bott_forc;

    // compute max level on block (maxval fortran function)
    int max_levels = 0;
    for (int jc = start_index; jc <= end_index; jc++)
        if (p_patch.dolic_c(blockNo, jc) > max_levels)
            max_levels = p_patch.dolic_c(blockNo, jc);

    // Initialize diagnostics and calculate mixing length scale
    for (int level = 0; level < p_constant.nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_cvmix.tke_Twin(blockNo, level, jc) = 0.0;
            p_internal.sqrttke(level, jc) = sqrt(max(0.0, p_internal.tke_old(level, jc)));
            p_cvmix.tke_Lmix(blockNo, level, jc) = sqrt(2.0) * p_internal.sqrttke(level, jc) /
                                    sqrt(max(1.0e-12, p_internal.Nsqr(level, jc)));
        }
    }

    if (p_constant_tke.tke_mxl_choice == 2) {
        calc_mxl_2<T>(blockNo, start_index, end_index, max_levels, p_constant_tke.mxl_min,
                           p_patch.dolic_c, p_cvmix.tke_Lmix, p_internal.dzw_stretched);
    } else if (p_constant_tke.tke_mxl_choice == 3) {
        // TODO(EnricoDeg): not default
    } else {
        // Error
    }

    // calculate diffusivities
    calc_diffusivity<T>(blockNo, start_index, end_index, max_levels, &p_constant_tke,
                     p_patch.dolic_c, p_cvmix.tke_Lmix, p_internal.sqrttke,
                     p_internal.Nsqr, p_internal.Ssqr,
                     p_internal.tke_Av, p_internal.tke_kv, p_cvmix.tke_Pr);

    // tke forcing
    calc_forcing<T>(blockNo, start_index, end_index, max_levels, p_constant.l_lc, p_constant_tke.only_tke,
                 p_patch.dolic_c, p_internal.Ssqr, p_internal.Nsqr, p_internal.tke_Av,
                 p_internal.tke_kv, p_cvmix.tke_Tspr, p_cvmix.tke_Tbpr,
                 p_cvmix.tke_plc, p_cvmix.tke_Tiwf, p_internal.forc);

    // vertical dissipation and diffusion solved implicitly
    build_diffusion_dissipation_tridiag<T>(blockNo, start_index, end_index, max_levels,
                                        p_patch.dolic_c, p_constant_tke.alpha_tke,
                                        p_internal.tke_Av, p_internal.dzt_stretched,
                                        p_internal.dzw_stretched,
                                        p_internal.ke, p_internal.a_dif, p_internal.b_dif,
                                        p_internal.c_dif);

    // copy tke_old
    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.tke_upd(level, jc) = p_internal.tke_old(level, jc);

    // upper boundary condition
    if (p_constant_tke.use_ubound_dirichlet) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_internal.sqrttke(0, jc) = 0.0;
            p_internal.forc(0, jc) = 0.0;
            tke_surf = max(p_constant_tke.tke_surf_min,
                           p_constant_tke.cd * p_internal.forc_tke_surf_2D(jc));
            p_internal.tke_upd(0, jc) = tke_surf;
            diff_surf_forc = p_internal.a_dif(1, jc) * tke_surf;
            p_internal.forc(1, jc) += diff_surf_forc;
            p_internal.a_dif(1, jc) = 0.0;
            p_internal.b_dif(0, jc) = 0.0;
            p_internal.c_dif(0, jc) = 0.0;
        }
    } else {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_internal.forc(0, jc) += (p_constant_tke.cd * pow(p_internal.forc_tke_surf_2D(jc), 1.5)) /
                                      p_internal.dzt_stretched(0, jc);
            p_internal.b_dif(0, jc) = p_internal.ke(0, jc) /
                                      (p_internal.dzt_stretched(0, jc) * p_internal.dzw_stretched(0, jc));
            diff_surf_forc = 0.0;
        }
    }

    // lower boundary condition
    if (p_constant_tke.use_lbound_dirichlet) {
        for (int jc = start_index; jc <= end_index; jc++) {
            int dolic = p_patch.dolic_c(blockNo, jc);
            p_internal.sqrttke(dolic, jc) = 0.0;
            p_internal.forc(dolic, jc) = 0.0;
            tke_bott = p_constant_tke.tke_min;
            p_internal.tke_upd(dolic, jc) = tke_bott;
            diff_bott_forc = p_internal.c_dif(dolic-1, jc) * tke_bott;
            p_internal.forc(dolic-1, jc) += diff_bott_forc;
            p_internal.c_dif(dolic-1, jc) = 0.0;
            p_internal.b_dif(dolic, jc) = 0.0;
            p_internal.a_dif(dolic, jc) = 0.0;
        }
    } else {
        for (int jc = start_index; jc <= end_index; jc++) {
            int dolic = p_patch.dolic_c(blockNo, jc);
            p_internal.b_dif(dolic, jc) = p_internal.ke(dolic-1, jc) /
                                            (p_internal.dzt_stretched(dolic, jc) *
                                            p_internal.dzw_stretched(dolic-1, jc));
            diff_bott_forc = 0.0;
        }
    }

    // construct tridiagonal matrix to solve diffusion and dissipation implicitely
    build_tridiag<T>(blockNo, start_index, end_index, max_levels, p_patch.dolic_c,
                  p_constant.dtime, p_constant_tke.c_eps, p_constant.nlevs,
                  p_internal.a_dif, p_internal.b_dif, p_internal.c_dif,
                  p_internal.sqrttke, p_cvmix.tke_Lmix, p_internal.tke_upd, p_internal.forc,
                  p_internal.a_tri, p_internal.b_tri, p_internal.c_tri, p_internal.d_tri);

    // solve the tri-diag matrix
    solve_tridiag<T>(blockNo, start_index, end_index, max_levels, p_patch.dolic_c,
                  p_internal.a_tri, p_internal.b_tri, p_internal.c_tri,
                  p_internal.d_tri, p_cvmix.tke, p_internal.cp, p_internal.dp);

    // diagnose implicite tendencies (only for diagnostics)
    // vertical diffusion of TKE
    tke_vertical_diffusion<T>(blockNo, start_index, end_index, max_levels, p_patch.dolic_c,
                           diff_surf_forc, diff_bott_forc,
                           p_internal.a_dif, p_internal.b_dif, p_internal.c_dif,
                           p_cvmix.tke, p_cvmix.tke_Tdif);

    // flux out of first box due to diffusion with Dirichlet boundary value of TKE
    // (tke_surf=tke_upd(1)) and TKE of box below (tke_new(2))
    if (p_constant_tke.use_ubound_dirichlet)
        tke_vertical_diffusion_ub_dirichlet<T>(blockNo, start_index, end_index, p_patch.dolic_c,
                                            tke_surf, p_internal.ke, p_internal.dzw_stretched,
                                            p_internal.dzt_stretched, p_cvmix.tke,
                                            p_cvmix.tke_Tdif);

    if (p_constant_tke.use_lbound_dirichlet)
        tke_vertical_diffusion_lb_dirichlet<T>(blockNo, start_index, end_index, p_patch.dolic_c,
                                            tke_bott, p_internal.ke, p_internal.dzw_stretched,
                                            p_internal.dzt_stretched, p_cvmix.tke,
                                            p_cvmix.tke_Tdif);

    // dissipation of TKE
    tke_vertical_dissipation<T>(blockNo, start_index, end_index, max_levels, p_patch.dolic_c,
                             p_constant.nlevs, p_constant_tke.c_eps, p_cvmix.tke_Lmix, p_internal.sqrttke,
                             p_cvmix.tke, p_cvmix.tke_Tdis);

    // reset tke to bounding values
    for (int level = 0; level < p_constant.nlevs+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.tke_unrest(level, jc) = p_cvmix.tke(blockNo, level, jc);

    // restrict values of TKE to tke_min, if IDEMIX is not used
    if (p_constant_tke.only_tke) {
        for (int level = 0; level < max_levels+1; level++)
            for (int jc = start_index; jc <= end_index; jc++)
                if (level < p_patch.dolic_c(blockNo, jc) + 1)
                    p_cvmix.tke(blockNo, level, jc) = max(p_cvmix.tke(blockNo, level, jc), p_constant_tke.tke_min);
    }

    // assign diagnostic variables
    for (int level = 0; level < max_levels+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < p_patch.dolic_c(blockNo, jc) + 1) {
                p_cvmix.tke_Tbpr(blockNo, level, jc) *= -1.0;
                p_cvmix.tke_Tbck(blockNo, level, jc) = (p_cvmix.tke(blockNo, level, jc) -
                                                       p_internal.tke_unrest(level, jc)) /
                                                       p_constant.dtime;
            }
        }
    }

    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_cvmix.tke_Tbck(blockNo, level, jc) = (p_cvmix.tke(blockNo, level, jc) -
                                                   p_internal.tke_unrest(level, jc)) /
                                                   p_constant.dtime;

    if (p_constant_tke.use_ubound_dirichlet) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_cvmix.tke_Twin(blockNo, 0, jc) = (p_cvmix.tke(blockNo, 0, jc) - p_internal.tke_old(0, jc)) /
                                               p_constant.dtime - p_cvmix.tke_Tdif(blockNo, 0, jc);
            p_cvmix.tke_Tbck(blockNo, 0, jc) = 0.0;
        }
    } else {
        for (int jc = start_index; jc <= end_index; jc++)
            p_cvmix.tke_Twin(blockNo, 0, jc) = (p_constant_tke.cd * pow(p_internal.forc_tke_surf_2D(jc), 1.5)) /
                                               p_internal.dzt_stretched(0, jc);
    }

    if (p_constant_tke.use_lbound_dirichlet) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (p_patch.dolic_c(blockNo, jc) > 0) {
                int dolic = p_patch.dolic_c(blockNo, jc);
                p_cvmix.tke_Twin(blockNo, dolic, jc) = (p_cvmix.tke(blockNo, dolic, jc) -
                                                       p_internal.tke_old(dolic, jc)) /
                                                       p_constant.dtime -
                                                       p_cvmix.tke_Tdif(blockNo, dolic, jc);
                p_cvmix.tke_Tbck(blockNo, dolic, jc) = 0.0;
            }
        }
    } else {
        for (int jc = start_index; jc <= end_index; jc++)
            if (p_patch.dolic_c(blockNo, jc) > 0)
                p_cvmix.tke_Twin(blockNo, p_patch.dolic_c(blockNo, jc), jc) = 0.0;
    }

    for (int level = 0; level < p_constant.nlevs+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_cvmix.tke_Ttot(blockNo, level, jc) = (p_cvmix.tke(blockNo, level, jc) -
                                                   p_internal.tke_old(level, jc)) / p_constant.dtime;

    for (int level = 0; level < p_constant.nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level >= p_patch.dolic_c(blockNo, jc)+1) {
                p_cvmix.tke_Lmix(blockNo, level, jc) = 0.0;
                p_cvmix.tke_Pr(blockNo, level, jc) = 0.0;
            }
        }
    }

    // the rest is for debugging
    for (int level = 0; level < p_constant.nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_cvmix.cvmix_dummy_1(blockNo, level, jc) = p_internal.tke_kv(level, jc);
            p_cvmix.cvmix_dummy_2(blockNo, level, jc) = p_internal.tke_Av(blockNo, level, jc);
            p_cvmix.cvmix_dummy_3(blockNo, level, jc) = p_internal.Nsqr(level, jc);
        }
    }
}

template <class T>
void calc_impl_cells(int blockNo, int start_index, int end_index,
                     t_patch_view<T, cpu_memview::mdspan, cpu_memview::dextents> p_patch,
                     t_cvmix_view<T, cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
                     t_ocean_state_view<cpu_memview::mdspan, cpu_memview::dextents> ocean_state,
                     t_atmo_fluxes_view<cpu_memview::mdspan, cpu_memview::dextents> atmos_fluxes,
                     t_atmos_for_ocean_view<T, cpu_memview::mdspan, cpu_memview::dextents> p_as,
                     t_sea_ice_view<T, cpu_memview::mdspan, cpu_memview::dextents> p_sea_ice,
                     t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
                     t_constant p_constant,
                     t_constant_tke p_constant_tke) {
    // initialization
    for (int level = 0; level < p_constant.nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_internal.tke_kv(level, jc) = 0.0;
            p_internal.tke_Av(blockNo, level, jc) = 0.0;
            if (p_constant.vert_mix_type == p_constant.vmix_idemix_tke) {
                p_cvmix.tke_Tiwf(blockNo, level, jc) = -1.0 * p_cvmix.iwe_Tdis(blockNo, level, jc);
            } else {
                p_cvmix.tke_Tiwf(blockNo, level, jc) = 0.0;
            }
            p_internal.dzt_stretched(level, jc) = p_patch.prism_center_dist_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);
        }
    }

    for (int level = 0; level < p_constant.nlevs; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.dzw_stretched(level, jc) = p_patch.prism_thick_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);

    // pre-integration
    for (int level = 0; level < p_constant.nlevs+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.tke_old(level, jc) = p_cvmix.tke(blockNo, level, jc);

    for (int jc = start_index; jc <= end_index; jc++) {
        T tau_abs = (1.0 - p_sea_ice.concsum(blockNo, jc))
                          * sqrt((atmos_fluxes.stress_xw(blockNo, jc) *
                                  atmos_fluxes.stress_xw(blockNo, jc))
                               + (atmos_fluxes.stress_yw(blockNo, jc) *
                                  atmos_fluxes.stress_yw(blockNo, jc)));
        p_internal.forc_tke_surf_2D(jc) = tau_abs / p_constant.OceanReferenceDensity;
    }

    // compute max level on block (maxval fortran function)
    int max_levels = 0;
    for (int jc = start_index; jc <= end_index; jc++)
        if (p_patch.dolic_c(blockNo, jc) > max_levels)
            max_levels = p_patch.dolic_c(blockNo, jc);

    // Loop over internal interfaces, surface (jk=1) and bottom (jk=kbot+1) excluded
    for (int level = 1; level < max_levels; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            if (level < p_patch.dolic_c(blockNo, jc)) {
                T rho_up = calculate_density(ocean_state.temp(blockNo, level-1, jc),
                                                  ocean_state.salt(blockNo, level-1, jc),
                                                  p_patch.zlev_i(level) *
                                                  p_constant.ReferencePressureIndbars);
                T rho_down = calculate_density(ocean_state.temp(blockNo, level, jc),
                                                    ocean_state.salt(blockNo, level, jc),
                                                    p_patch.zlev_i(level) *
                                                    p_constant.ReferencePressureIndbars);
                p_internal.Nsqr(level, jc) = p_constant.grav / p_constant.OceanReferenceDensity *
                                             (rho_down - rho_up) *
                                             p_patch.inv_prism_center_dist_c(blockNo, level, jc) /
                                             ocean_state.stretch_c(blockNo, jc);
                p_internal.Ssqr(level, jc) = pow((ocean_state.p_vn_x1(blockNo, level-1, jc) -
                                             ocean_state.p_vn_x1(blockNo, level, jc) ) *
                                             p_patch.inv_prism_center_dist_c(blockNo, level, jc) /
                                             ocean_state.stretch_c(blockNo, jc) , 2.0) +
                                             pow((ocean_state.p_vn_x2(blockNo, level-1, jc) -
                                             ocean_state.p_vn_x2(blockNo, level, jc) ) *
                                             p_patch.inv_prism_center_dist_c(blockNo, level, jc) /
                                             ocean_state.stretch_c(blockNo, jc) , 2.0) +
                                             pow((ocean_state.p_vn_x3(blockNo, level-1, jc) -
                                             ocean_state.p_vn_x3(blockNo, level, jc) ) *
                                             p_patch.inv_prism_center_dist_c(blockNo, level, jc) /
                                             ocean_state.stretch_c(blockNo, jc) , 2.0);
            }
        }
    }

    for (int level = 0; level < max_levels; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.dzw_stretched(level, jc) = p_patch.prism_thick_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);
    for (int level = 0; level < max_levels+1; level++)
        for (int jc = start_index; jc <= end_index; jc++)
            p_internal.dzt_stretched(level, jc) = p_patch.prism_center_dist_c(blockNo, level, jc) *
                                                  ocean_state.stretch_c(blockNo, jc);

    // integration
    integrate<T>(blockNo, start_index, end_index, p_patch, p_cvmix, p_internal, p_constant,
              p_constant_tke);

    //  write tke vert. diffusivity to vert tracer diffusivities
    for (int level = 0; level < p_constant.nlevs+1; level++) {
        for (int jc = start_index; jc <= end_index; jc++) {
            p_cvmix.a_temp_v(blockNo, level, jc) = p_internal.tke_kv(level, jc);
            p_cvmix.a_salt_v(blockNo, level, jc) = p_internal.tke_kv(level, jc);
        }
    }
}

#endif  // SRC_BACKENDS_CPU_CPU_KERNELS_HPP_
