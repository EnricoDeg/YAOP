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

void TKE::calc(int start_block, int end_block,
               double *depth_CellInterface, double *prism_center_dist_c,
               double *inv_prism_center_dist_c, double *prism_thick_c, 
               int *dolic_c, int *dolic_e, double *zlev_i, double *wet_c,
               int *edges_cell_idx, int *edges_cell_blk,
               double *temp, double *salt, double *stretch_c, double *eta_c,
               double *tke, double *tke_plc_in, double *hlc_in, double *wlc_in,
               double *u_stokes_in, double *a_veloc_v, double *a_temp_v, double *a_salt_v,
               double *iwe_Tdis, double *cvmix_dummy_1, double *cvmix_dummy_2,
               double *cvmix_dummy_3, double *tke_Tbpr, double *tke_Tspr,
               double *tke_Tdif, double *tke_Tdis, double *tke_Twin,
               double *tke_Tiwf, double *tke_Tbck, double *tke_Ttot,
               double *tke_Lmix, double *tke_Pr, double *stress_xw,
               double *stress_yw, double *fu10, double *concsum) {

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
