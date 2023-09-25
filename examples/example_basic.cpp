#include <memory>
#include "TKE.hpp"

int main(int argc, char ** argv) {
    std::shared_ptr<TKE> ocean_physics;
    ocean_physics.reset(new TKE());

    ocean_physics.reset();
}