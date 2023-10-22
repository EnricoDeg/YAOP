#include <memory>
extern "C" {
#include "TKE.h"
}
#include "TKE.hpp"

static std::unique_ptr<TKE> impl = nullptr;

void TKE_Init(int nproma, int nlevs, int nblocks,
              int block_size, int start_index, int end_index) {
    impl.reset(new TKE(nproma, nlevs, nblocks,
                       block_size, start_index, end_index));
}

void TKE_Finalize() {
    impl.reset();
}

void TKE_Calc(int start_block, int end_block, double * temperature) {
    impl->calc(start_block, end_block, temperature);
}
