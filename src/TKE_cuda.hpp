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

#ifndef TKE_CUDA_HPP
#define TKE_CUDA_HPP

#include "src/TKE_backend.hpp"

class TKE_cuda : public TKE_backend {

    public:
        TKE_cuda(int nproma, int nlevs, int nblocks,
                 int block_size, int start_index, int end_index);
        ~TKE_cuda();

    protected:
        void calc_impl(int start_block, int end_block, struct t_patch p_patch, struct t_cvmix p_cvmix);

    private:
        bool is_view_init;

};

#endif /* TKE_CUDA_HPP */
