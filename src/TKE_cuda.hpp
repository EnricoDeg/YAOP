#ifndef TKE_CUDA_HPP
#define TKE_CUDA_HPP

#include "TKE_backend.hpp"

class TKE_cuda : public TKE_backend {

    public:
        TKE_cuda(int nproma, int nlevs, int nblocks,
                 int block_size, int start_index, int end_index);
        ~TKE_cuda();

    protected:
        void calc_impl(int start_block, int end_block, double *tke);

    private:
        bool is_view_init;

};

#endif /* TKE_CUDA_HPP */
