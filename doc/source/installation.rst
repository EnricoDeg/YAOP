.. _installation:

Installation
============

Clone the repository::

  git clone https://github.com/EnricoDeg/TKE.git

Init and update submodules::

  git submodule init
  git submodule update

The Cmake configuration has the following options (default is disabled):

 - ENABLE_C: enables the C bindings and the examples and tests based on them

 - ENABLE_FORTRAN: enables the Fortran bindings and the examples and tests based on them

 - ENABLE_CUDA: enable the CUDA backend of the GPU implementation (CPU implementation not compiled)

 - ENABLE_HIP: enable the HIP backend of the GPU implementation (CPU implementation not compiled)

 - ENABLE_EXAMPLES: compile files in ``examples`` folder

 - ENABLE_TESTS: install gtest and compile files in ``tests`` folder

The default installation is straightforeward::

  mkdir build
  cd build
  cmake ..
  make
  make install

And it is compiling only the files in the ``src`` folder for the CPU implementation.

