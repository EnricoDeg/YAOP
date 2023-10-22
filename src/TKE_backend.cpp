#include <iostream>
#include "TKE_backend.hpp"

TKE_backend::TKE_backend(int nproma, int nlevs, int nblocks,
                         int block_size, int start_index, int end_index) 
    : m_nproma(nproma), m_nlevs(nlevs), m_nblocks(nblocks),
      m_block_size(block_size), m_start_index(start_index), m_end_index(end_index) {
}

void TKE_backend::calc(int start_block, int end_block, double *tke) {
    this->calc_impl(start_block, end_block, tke);
}
