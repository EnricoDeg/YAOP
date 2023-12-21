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

#include "src/shared/interface/data_struct.hpp"

void fill_struct(struct t_atmo_fluxes *atmo_fluxes, double *stress_xw, double *stress_yw) {
    atmo_fluxes->stress_xw = stress_xw;
    atmo_fluxes->stress_yw = stress_yw;
}

void fill_struct(struct t_atmos_for_ocean *p_as, double *fu10) {
    p_as->fu10 = fu10;
}

void fill_struct(struct t_sea_ice *p_sea_ice, double *concsum) {
    p_sea_ice->concsum = concsum;
}
