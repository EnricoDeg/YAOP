#include  <experimental/mdspan>
#include <iostream>


namespace stdex = std::experimental;

template<class T>
class TKE {
    public:
        
    	TKE(int nproma, int nlevs, int nblocks);
        void calc(T * temperature);

    protected:
    private:
        int m_nproma;
        int m_nlevs;
        int m_nblocks;
        stdex::mdspan<T, stdex::dextents<size_t, 3>> m_temperature;


        template <
            class ExtsA, class LayA, class AccA
        >
        void print_field(
        stdex::mdspan<T, ExtsA, LayA, AccA> a
        ); // requires ExtsA::rank() == 3
};
