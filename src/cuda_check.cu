#include <iostream>
#include "cuda_check.hpp"

void check(cudaError_t err) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
        std::exit(-1);
    }
}
