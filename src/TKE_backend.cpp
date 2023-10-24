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

#include "src/TKE_backend.hpp"
#include <iostream>

TKE_backend::TKE_backend(int nproma, int nlevs, int nblocks,
                         int block_size, int start_index, int end_index)
    : m_nproma(nproma), m_nlevs(nlevs), m_nblocks(nblocks),
      m_block_size(block_size), m_start_index(start_index), m_end_index(end_index) {
}

void TKE_backend::calc(int start_block, int end_block, struct t_patch p_patch, struct t_cvmix p_cvmix) {
    this->calc_impl(start_block, end_block, p_patch, p_cvmix);
}
