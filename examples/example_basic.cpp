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

#include <memory>
#include "src/TKE.hpp"

int main(int argc, char ** argv) {
    int nproma = 1200;
    int nlevs = 64;
    int nblocks = 2;
    int vert_mix_type = 2;
    int vmix_idemix_tke = 4;
    int vert_cor_type = 0;
    double dtime = 0.0;
    double OceanReferenceDensity = 1025.022;
    double grav = 9.80665;
    int l_lc = 0;
    double clc = 0.15;
    double ReferencePressureIndbars = 1035.0*grav*1.0e-4;
    double pi = 3.14159265358979323846264338327950288;

    std::shared_ptr<TKE> ocean_physics;
    ocean_physics.reset(new TKE(nproma, nlevs, nblocks, vert_mix_type, vmix_idemix_tke,
                                vert_cor_type, dtime, OceanReferenceDensity, grav,
                                l_lc, clc, ReferencePressureIndbars, pi));

    ocean_physics.reset();
}
