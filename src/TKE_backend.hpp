#ifndef TKE_BACKEND_HPP
#define TKE_BACKEND_HPP

#include <memory>

class TKE_backend {

    public:
        typedef std::shared_ptr<TKE_backend> Ptr;
        TKE_backend(int nproma, int nlevs, int nblocks);
        virtual ~TKE_backend() = default;
        void calc();

    protected:
        virtual void calc_impl() = 0;

    protected:
        double *rho_up;
        double *rho_down;
        int m_nproma;
        int m_nlevs;
        int m_nblocks;

};

#endif /* TKE_BACKEND_HPP */
