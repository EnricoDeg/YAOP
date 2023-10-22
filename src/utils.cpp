#include "utils.hpp"

void get_index_range(int subset_block_size, int subset_start_block, int subset_end_block,
                     int subset_start_index, int subset_end_index, int current_block, 
                     int *start_index, int *end_index) {
    *start_index = 0;
    *end_index = subset_block_size;

    if (current_block == subset_start_block)
        *start_index = subset_start_index;
    if (current_block == subset_end_block)
        *end_index = subset_end_index;
}

