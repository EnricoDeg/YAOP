#include <iostream>
#include  <experimental/mdspan>
#include "TKE.hpp"
#include "TKE_backend.hpp"
#include "TKE_cuda.hpp"

namespace stdex = std::experimental;

struct TKE::Impl {
  TKE_backend::Ptr backend;
};

TKE::TKE(int nproma, int nlevs, int nblocks, 
         int block_size, int start_index, int end_index) 
    : m_impl(new Impl) {
    m_nproma = nproma;
    m_nlevs = nlevs;
    m_nblocks = nblocks;
    std::cout << "Setting nproma to " << m_nproma << std::endl;
    std::cout << "Setting nlevs to " << m_nlevs << std::endl;
    std::cout << "Setting nblocks to " << m_nblocks << std::endl;

    m_impl->backend = TKE_backend::Ptr(new TKE_cuda(nproma, nlevs, nblocks, 
                                                    block_size, start_index, end_index));

}

TKE::~TKE() {
    std::cout << "Finalizing TKE... " << std::endl;
    delete m_impl;
}

void TKE::calc(int start_block, int end_block, double *tke, int *dolic_c) {
    m_impl->backend->calc(start_block, end_block, tke, dolic_c);
}

template <
    class T,
    class ExtsA, class LayA, class AccA
>
void TKE::print_field(
    stdex::mdspan<T, ExtsA, LayA, AccA> a
    ) // requires ExtsA::rank() == 3
{
    for(int i = 0; i < a.extent(0); ++i)
        for(int j = 0; j < a.extent(1); ++j)
            for(int k = 0; k < a.extent(2); ++k)
                std::cout << "field(" << i << "," << j << "," << k << ") = " << a(i,j,k) << std::endl;
}
