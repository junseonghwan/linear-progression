# Linear Progression Model

[![Build Status](https://travis-ci.com/junseonghwan/linear-progression.svg?token=wxZvzvzdwz1aU7zpr7vw&branch=master)](https://travis-ci.com/junseonghwan/linear-progression)

The core library is written in C++. The library can be called from Python and the sample code is provided (see `run_xxx.py`). The python code depends on `ctypes` to load the library found in `bin/liblpm_lib.dylib`.

## To compile from the source
Requirements:
+ `cmake`, `gsl`, `boost`, `libomp`
+ Install [SPF](https://github.com/junseonghwan/spf)
+ `git clone https://github.com/junseonghwan/linear-progression.git`
+ `mkdir linear-progression/build`
+ `cd linear-progression/build`
+ `cmake ..`
+ `make install`

To develop the code
+ `mkdir linear-progression/xcode`
+ `cd linear-progression/xcode`
+ `cmake .. -GXcode`
+ This wil create `lpm.xcodeproj`, which can be loaded by Xcode.

TODO:
+ Provide unified interface. Users should only need to provide a config file and mode of execution (i.e., one of model selection, PG). Python program should read in the config file and perform the tasks similar to `run_xxx.py`.
