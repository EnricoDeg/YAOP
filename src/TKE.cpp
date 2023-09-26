#include <iostream>
#include  <experimental/mdspan>
#include "TKE.hpp"

namespace stdex = std::experimental;

template <class T>
TKE<T>::TKE(int nproma, int nlevs, int nblocks) {
    m_nproma = nproma;
    m_nlevs = nlevs;
    m_nblocks = nblocks;
    std::cout << "Setting nproma to " << m_nproma << std::endl;
    std::cout << "Setting nlevs to " << m_nlevs << std::endl;
    std::cout << "Setting nblocks to " << m_nblocks << std::endl;
}

template <class T>
void TKE<T>::calc(T *temperature) {
    m_temperature = stdex::mdspan<T, stdex::dextents<size_t, 3>>(temperature, stdex::extents{m_nblocks, m_nlevs, m_nproma});

    print_field(m_temperature);
}

template <class T>
template <
    class ExtsA, class LayA, class AccA
>
void TKE<T>::print_field(
    stdex::mdspan<T, ExtsA, LayA, AccA> a
    ) // requires ExtsA::rank() == 3
{
    for(int i = 0; i < a.extent(0); ++i)
        for(int j = 0; j < a.extent(1); ++j)
            for(int k = 0; k < a.extent(2); ++k)
                std::cout << "field(" << i << "," << j << "," << k << ") = " << a(i,j,k) << std::endl;
}

template class TKE<float>;
template class TKE<double>;
