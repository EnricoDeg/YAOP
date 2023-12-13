.. _library:

Library
=======

The library is designed to have three steps:

 - Initialization: this happens before the time loop and it provides the main information about the model which are needed inside the library. Internally all the arrays are allocated.

 - Computation: this happens inside the time loop. The provided pointers are not expected to change during the time loop.

 - Finalization: this happens after the time loop. The internal memory is cleaned.

The library is written in C++ but C and Fortran interfaces are provided. The `examples` folder provides simple benchmarks of the vertical mixing scheme running standalone and with different drivers written in the supported languages.

.. doxygenfile:: TKE.hpp
   :project: TKE

.. toctree::
   :maxdepth: 2

   bindings
