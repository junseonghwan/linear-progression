import numpy as np 
import ctypes
import os

lpm_lib = np.ctypeslib.load_library("liblpm_lib.dylib", "build/src/")

curr_dir = os.getcwd()
data_path = curr_dir + "/data/raphael/25genes/error0.01/rep0/"
input_path = data_path + "matrix.csv"

seed = ctypes.c_long(1)
_input_path = ctypes.create_string_buffer(input_path.encode())
_data_path = ctypes.create_string_buffer(data_path.encode())
model_len = ctypes.c_uint(5)
n_pg_iter = ctypes.c_uint(10)
n_particles = ctypes.c_uint(100)
n_smc_iter = ctypes.c_uint(100)
n_kernel_iter = ctypes.c_uint(1)
n_mh_w_gibbs_iter = ctypes.c_uint(20)
has_passenger = ctypes.c_bool(False)
swap_prob = ctypes.c_double(0.2)
fbp_max = ctypes.c_double(0.0)
bgp_max = ctypes.c_double(0.3)

lpm_lib.run_pg(seed, _input_path, _data_path, model_len, n_pg_iter, 
                n_particles, n_smc_iter, n_kernel_iter, n_mh_w_gibbs_iter,
                has_passenger, swap_prob, fbp_max, bgp_max)

