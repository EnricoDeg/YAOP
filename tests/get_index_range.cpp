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
#include "src/utils.hpp"

// Check the default value assignment for get_index_range
TEST(get_index_range, default_value) {
  int start_index, end_index;
  start_index = -1;
  end_index = -1;

  // Assert that the default start and end index values are assigned
  get_index_range(1, -1, -1, 0, 2, 0, &start_index, &end_index);
  ASSERT_EQ(start_index, 0);
  ASSERT_EQ(end_index, 1);
}

// Check the default value assignment for get_index_range
TEST(get_index_range, start_index_assignment) {
  int start_index, end_index;
  start_index = -1;
  end_index = -1;

  // Assert that the default start and end index values are assigned
  get_index_range(-1, 1, -1, 2, 2, 1, &start_index, &end_index);
  ASSERT_EQ(start_index, 2);
  ASSERT_EQ(end_index, -1);
}

// Check the default value assignment for get_index_range
TEST(get_index_range, end_index_assignment) {
  int start_index, end_index;
  start_index = -1;
  end_index = -1;

  // Assert that the default start and end index values are assigned
  get_index_range(1, -1, 2, 0, 2, 2, &start_index, &end_index);
  ASSERT_EQ(start_index, 0);
  ASSERT_EQ(end_index, 2);
}
