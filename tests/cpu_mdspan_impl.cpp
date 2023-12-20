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

#include <gtest/gtest.h>
#include "src/backends/CPU/cpu_memory.hpp"

// Test memory view creating from allocated memory
TEST(cpu_mdspan_impl, memview) {
    int nblocks = 1;
    int nlevs = 1;
    int nproma = 1;

    double *test_ptr = reinterpret_cast<double *>(malloc(nblocks * nlevs * nproma * sizeof(double)));
    mdspan_3d<double> test = cpu_mdspan_impl<double>::memview(test_ptr, nblocks, nlevs, nproma);
    ASSERT_EQ(test.size(), nblocks*nlevs*nproma);
    free(test_ptr);

    test_ptr = NULL;
    mdspan_3d<double> test_null = cpu_mdspan_impl<double>::memview(test_ptr, nblocks, nlevs, nproma);
    ASSERT_EQ(test_null.size(), 0);
}
