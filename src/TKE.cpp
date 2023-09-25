#include <iostream>
#include "TKE.hpp"

TKE::TKE(int nproma, int nlevs, int nblocks) {
	m_nproma = nproma;
	m_nlevs = nlevs;
	m_nblocks = nblocks;
    std::cout << "Setting nproma to " << m_nproma << std::endl;
    std::cout << "Setting nlevs to " << m_nlevs << std::endl;
    std::cout << "Setting nblocks to " << m_nblocks << std::endl;
}

