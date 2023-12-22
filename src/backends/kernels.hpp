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

#include "src/shared/constants/constants_thermodyn.hpp"

#if defined CUDA || defined HIP

#else
#include <cmath>
#include <algorithm>
using std::max;
#endif

/*! \brief Compute pointwise density from temperature, salinity and pressure.
*
*/
#if defined CUDA || defined HIP
__device__
#endif
template <class T>
inline
T calculate_density(T temp, T salt, T pressure) {
    T dvs, fne, fst, qn3;
    T qnq, qvs, s, s3h;
    T t, denom, s__2;
    T rho;
    T zero = 0.0;
    T two = 2.0;

    // This is the adisit part, that transforms potential in in-situ temperature
    qnq = -pressure * (-a_a3 + pressure * a_c3);
    qn3 = -pressure * a_a4;
    qvs = (pressure * (a_b1 - a_d * pressure)) *
          (salt - z_sref) + pressure * (a_a1 + pressure * (a_c1 - a_e1 * pressure));

    dvs = (a_b2 * pressure) * (salt - z_sref) +
           1.0 + pressure * (-a_a2 + pressure * (a_c2 - a_e2 * pressure));

    t   = (temp + qvs) / dvs;
    fne = - qvs + t * (dvs + t * (qnq + t * qn3)) - temp;

    fst = dvs + t * (2.0 * qnq + 3.0 * qn3 * t);

    t    = t - fne / fst;
    s    = max(salt, zero);
    s__2 = pow(s, two);
    s3h  = s * sqrt(s);

    rho = r_a0 + t * (r_a1 + t * (r_a2 + t * (r_a3 + t * (r_a4 + t * r_a5))))
        + s * (r_b0 + t * (r_b1 + t * (r_b2 + t * (r_b3 + t * r_b4))))
        + r_d0 * s__2 + s3h * (r_c0 + t * (r_c1 + r_c2 * t));

    denom = 1.0 - pressure / (pressure * (r_h0 + t *
            (r_h1 + t * (r_h2 + t * r_h3))
            + s * (r_ai0 + t * (r_ai1 + r_ai2 * t))
            + r_aj0 * s3h + (r_ak0 + t * (r_ak1 + t * r_ak2)
            + s * (r_am0 + t * (r_am1 + t * r_am2))) * pressure)
            + r_e0 + t * (r_e1 + t * (r_e2 + t * (r_e3 + t * r_e4)))
            + s * (r_f0 + t * (r_f1 + t * (r_f2 + t * r_f3)))
            + s3h * (r_g0 + t * (r_g1 + r_g2 * t)));

    return rho/denom;
}

#endif  // SRC_BACKENDS_KERNELS_HPP_
