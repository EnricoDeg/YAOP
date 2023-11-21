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

#include "src/utils.hpp"

void get_index_range(int subset_block_size, int subset_start_block, int subset_end_block,
                     int subset_start_index, int subset_end_index, int current_block,
                     int *start_index, int *end_index) {
    *start_index = 0;
    *end_index = subset_block_size - 1;

    if (current_block == subset_start_block)
        *start_index = subset_start_index;
    if (current_block == subset_end_block)
        *end_index = subset_end_index;
}
