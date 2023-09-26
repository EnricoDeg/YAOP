#include  <experimental/mdspan>
#include <iostream>


namespace stdex = std::experimental;

class TKE {
    public:
        
    	TKE(int nproma, int nlevs, int nblocks);
        void calc(double * temperature);

    protected:
    private:
        int m_nproma;
        int m_nlevs;
        int m_nblocks;
        stdex::mdspan<double, stdex::dextents<size_t, 3>> m_temperature;


        template <
            class T,
            class ExtsA, class LayA, class AccA
        >
        void print_field(
        stdex::mdspan<T, ExtsA, LayA, AccA> a
        ); // requires ExtsA::rank() == 3
};
