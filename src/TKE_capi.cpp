#include <memory>
extern "C" {
#include "TKE.h"
}
#include "TKE.hpp"

static std::unique_ptr<TKE> impl = nullptr;

void TKE_Init() {
    impl.reset(new TKE());
}