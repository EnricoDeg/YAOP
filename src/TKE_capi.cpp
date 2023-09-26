#include <memory>
extern "C" {
#include "TKE.h"
}
#include "TKE.hpp"

static std::unique_ptr<TKE<double>> impl = nullptr;

void TKE_Init(int nproma, int nlevs, int nblocks) {
    impl.reset(new TKE<double>(nproma, nlevs, nblocks));
}

void TKE_Calc(double * temperature) {
    impl->calc(temperature);
}
