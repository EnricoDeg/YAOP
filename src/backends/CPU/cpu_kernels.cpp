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

#include <algorithm>
#include <cmath>
#include "src/backends/CPU/cpu_kernels.hpp"
#include "src/backends/kernels.hpp"
#include "src/shared/constants/constants_thermodyn.hpp"

void calc_impl_cells(int blockNo, int start_index, int end_index,
                     t_patch_view<cpu_memview::mdspan, cpu_memview::dextents> p_patch,
                     t_cvmix_view<cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
                     t_ocean_state_view<cpu_memview::mdspan, cpu_memview::dextents> ocean_state,
                     t_atmo_fluxes_view<cpu_memview::mdspan, cpu_memview::dextents> atmos_fluxes,
                     t_atmos_for_ocean_view<cpu_memview::mdspan, cpu_memview::dextents> p_as,
                     t_sea_ice_view<cpu_memview::mdspan, cpu_memview::dextents> p_sea_ice,
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
        double tau_abs = (1.0 - p_sea_ice.concsum(blockNo, jc))
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
                double rho_up = calculate_density(ocean_state.temp(blockNo, level-1, jc),
                                                  ocean_state.salt(blockNo, level-1, jc),
                                                  p_patch.zlev_i(level) *
                                                  p_constant.ReferencePressureIndbars);
                double rho_down = calculate_density(ocean_state.temp(blockNo, level, jc),
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
}

void integrate(int jc, int blockNo,
               t_patch_view<cpu_memview::mdspan, cpu_memview::dextents> p_patch,
               t_cvmix_view<cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
               t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
               t_constant p_constant,
               t_constant_tke p_constant_tke) {
}

void calc_impl_edges(int blockNo, int start_index, int end_index,
                     t_patch_view<cpu_memview::mdspan, cpu_memview::dextents> p_patch,
                     t_cvmix_view<cpu_memview::mdspan, cpu_memview::dextents> p_cvmix,
                     t_tke_internal_view<cpu_memview::mdspan, cpu_memview::dextents> p_internal,
                     t_constant p_constant) {
}
