#include <memory>
extern "C" {
#include "TKE.h"
}
#include "TKE.hpp"

static std::unique_ptr<TKE> impl = nullptr;

void TKE_Init(int nproma, int nlevs, int nblocks) {
    impl.reset(new TKE(nproma, nlevs, nblocks));
}

void TKE_Finalize() {
    impl.reset();
}

void TKE_Calc(double * temperature) {
    impl->calc(temperature);
}
