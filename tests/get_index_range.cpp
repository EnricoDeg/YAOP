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
#include "src/shared/utils.hpp"

// Test the default value assignment for start_index and end_index
TEST(get_index_range, default_value) {
    int start_index = -1;
    int end_index = -1;
    int current_block = 1;

    // default value that should be assigned to end_index
    int subset_block_size = 1;

    // These should be different than current_block
    int subset_start_block = -1;
    int subset_end_block = -1;

    get_index_range(subset_block_size, subset_start_block, subset_end_block,
                    -1, -1, current_block, &start_index, &end_index);

    // Assert that the default start and end index values are assigned
    ASSERT_EQ(start_index, 0);
    ASSERT_EQ(end_index, subset_block_size-1);
}

// Test that start_index is correctly assigned to subset_start_index if current_block equals subset_start_block
TEST(get_index_range, start_index_assignment) {
    int start_index = -1;
    int end_index = -1;

    // default value that should be assigned to end_index
    int subset_block_size = 1;

    // value that should be assigned to start_index
    int subset_start_index = 2;

    // these two values should be equal for the comparison inside get_index_range
    int current_block = 1;
    int subset_start_block = 1;

    // This should be different than current_block
    int subset_end_block = -1;

    get_index_range(subset_block_size, subset_start_block, subset_end_block,
                    subset_start_index, -1, current_block, &start_index, &end_index);

    // Assert that the start_index was assigned subset_start_index and end_index was assigned the default value
    ASSERT_EQ(start_index, subset_start_index);
    ASSERT_EQ(end_index, subset_block_size-1);
}

// Test that end_index is correctly assigned to subset_end_index if current_block equals subset_end_block
TEST(get_index_range, end_index_assignment) {
    int start_index = -1;
    int end_index = -1;

    // default value for end_index, should not be assined in this test
    int subset_block_size = 1;

    // value that should be assigned to end_index
    int subset_end_index = 2;

    // these two values should be equal for the comparison inside get_index_range
    int current_block = 1;
    int subset_end_block = 1;

    // This should be different than current_block
    int subset_start_block = -1;

    get_index_range(subset_block_size, subset_start_block, subset_end_block,
                    -1, subset_end_index, current_block, &start_index, &end_index);

    // Assert that the start_index was assigned the default value and end_index was assigned subset_end_index
    ASSERT_EQ(start_index, 0);
    ASSERT_EQ(end_index, subset_end_index);
}
