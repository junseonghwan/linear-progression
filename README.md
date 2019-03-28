# Linear Progression Model

The core library is written in C++. The library can be called from Python and the sample code is provided (see `run.py`). The python sample code uses `ctype` to load the library found in `bin/liblpm_lib.dylib`.

TODO:
+ Users should edit a config file. Python program should read in the config file and performs the tasks similar to `run.py` to call C++ library.

## To compile from the source
Requirements:
+ `cmake`, `gsl`, `boost`, `libomp`.
+ Install [SPF](https://github.com/junseonghwan/spf)