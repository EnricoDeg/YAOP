#ifndef TKE_HPP
#define TKE_HPP

#include <experimental/mdspan>
#include <iostream>


namespace stdex = std::experimental;

class TKE {
    public:
        
    	TKE(int nproma, int nlevs, int nblocks,
            int block_size, int start_index, int end_index);
        ~TKE();
        void calc(int start_block, int end_block, double *tke);

    protected:
    private:
        int m_nproma;
        int m_nlevs;
        int m_nblocks;
        struct Impl;
        Impl *m_impl;


        template <
            class T,
            class ExtsA, class LayA, class AccA
        >
        void print_field(
        stdex::mdspan<T, ExtsA, LayA, AccA> a
        ); // requires ExtsA::rank() == 3
};

#endif /* TKE_HPP */
