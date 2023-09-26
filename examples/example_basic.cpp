#include <memory>
#include "TKE.hpp"

int main(int argc, char ** argv) {
	int nproma = 1200;
	int nlevs = 64;
	int nblocks = 2;
    std::shared_ptr<TKE<double>> ocean_physics;
    ocean_physics.reset(new TKE<double>(nproma, nlevs, nblocks));

    ocean_physics.reset();
}
