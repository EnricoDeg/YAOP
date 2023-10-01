#include <iostream>
#include "TKE_backend.hpp"

TKE_backend::TKE_backend(int nproma, int nlevs, int nblocks) 
    : m_nproma(nproma), m_nlevs(nlevs), m_nblocks(nblocks) {
}

void TKE_backend::calc() {
    this->calc_impl();
}
