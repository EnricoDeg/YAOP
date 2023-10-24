# TKE

### Installation
`mkdir build`

`cd build`

`CC=nvc CXX=nvc++ FC=nvfortran CUDACXX=nvcc cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ..`

`make`

`make install`

### Development
Setup to use pre-commit:

`pip install pre-commit`

`pre-commit install`

`pip install cpplint`
