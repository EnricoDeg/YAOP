#ifndef TKE_CUDA_HPP
#define TKE_CUDA_HPP

#include "TKE_backend.hpp"

class TKE_cuda : public TKE_backend {

    public:
        TKE_cuda(int nproma, int nlevs, int nblocks);
        ~TKE_cuda();

    protected:
        void calc_impl();

};

#endif /* TKE_CUDA_HPP */
