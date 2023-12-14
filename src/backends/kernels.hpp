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

#ifndef SRC_BACKENDS_KERNELS_HPP_
#define SRC_BACKENDS_KERNELS_HPP_

/*! \brief Compute pointwise density from temperature, salinity and pressure.
*
*/
#if defined CUDA || defined HIP
__device__
#endif
double calculate_density(double temp, double salt, double pressure);

#endif  // SRC_BACKENDS_KERNELS_HPP_
