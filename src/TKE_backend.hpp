#ifndef TKE_BACKEND_HPP
#define TKE_BACKEND_HPP

#include <memory>
#include "data_struct.hpp"

class TKE_backend {

    public:
        typedef std::shared_ptr<TKE_backend> Ptr;
        TKE_backend(int nproma, int nlevs, int nblocks,
                    int block_size, int start_index, int end_index);
        virtual ~TKE_backend() = default;
        void calc(int start_block, int end_block, struct t_patch p_patch, struct t_cvmix p_cvmix);

    protected:
        virtual void calc_impl(int start_block, int end_block, struct t_patch p_patch, struct t_cvmix p_cvmix) = 0;

    protected:
        int m_nproma;
        int m_nlevs;
        int m_nblocks;
        int m_block_size;
        int m_start_index;
        int m_end_index;

        double *m_rho_up;
        double *m_rho_down;
        double *m_tke_old;
        double *m_tke_Av;
        double *m_tke_kv;
        double *m_tke_iw_alpha_c;
        double *m_tke_iwe;
        double *m_tke_iwe_forcing;
        double *m_forc_tke_surf_2D;
        double *m_forc_rho_surf_2D;
        double *m_bottom_fric_2D;
        double *m_s_c;
        double *m_dzw_stretched;
        double *m_dzt_stretched;
        double *m_pressure;
        double *m_Nsqr;
        double *m_Ssqr;

};

#endif /* TKE_BACKEND_HPP */
