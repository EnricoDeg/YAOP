#ifndef TKE_BACKEND_HPP
#define TKE_BACKEND_HPP

#include <memory>

class TKE_backend {

    public:
        typedef std::shared_ptr<TKE_backend> Ptr;
        TKE_backend(int nproma, int nlevs, int nblocks,
                    int block_size, int start_index, int end_index);
        virtual ~TKE_backend() = default;
        void calc(int start_block, int end_block, double *tke, int *dolic_c);

    protected:
        virtual void calc_impl(int start_block, int end_block, double *tke, int *dolic_c) = 0;

    protected:
        double *rho_up;
        double *rho_down;
        double *tke_old;
        int m_nproma;
        int m_nlevs;
        int m_nblocks;
        int m_block_size;
        int m_start_index;
        int m_end_index;

};

#endif /* TKE_BACKEND_HPP */
