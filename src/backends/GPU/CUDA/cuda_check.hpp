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

#ifndef SRC_BACKENDS_GPU_CUDA_CUDA_CHECK_HPP_
#define SRC_BACKENDS_GPU_CUDA_CUDA_CHECK_HPP_

#include <cuda_runtime.h>

/*! \brief Check error returned from CUDA function.
 *
 *  If not cudaSuccess, print a string describing the error and exit.
 */
void check(cudaError_t err);

#endif  // SRC_BACKENDS_GPU_CUDA_CUDA_CHECK_HPP_
