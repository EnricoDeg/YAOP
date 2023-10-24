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

#ifndef SRC_TKE_CUDA_HPP_
#define SRC_TKE_CUDA_HPP_

#include "src/TKE_backend.hpp"

class TKE_cuda : public TKE_backend {
 public:
    TKE_cuda(int nproma, int nlevs, int nblocks,
             int block_size, int start_index, int end_index);
    ~TKE_cuda();

 protected:
    void calc_impl(int start_block, int end_block, struct t_patch p_patch, struct t_cvmix p_cvmix,
                       struct t_ocean_state ocean_state, struct t_atmo_fluxes atmos_fluxes,
                       struct t_atmos_for_ocean p_as, struct t_sea_ice p_sea_ice);

 private:
    bool is_view_init;
};

#endif  // SRC_TKE_CUDA_HPP_
