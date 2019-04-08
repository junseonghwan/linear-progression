import numpy as np 
from scipy.special import logsumexp
import ctypes
import os
import platform

if platform.system() == "Linux":
    lpm_lib = np.ctypeslib.load_library("liblpm_lib.so", "bin/")    
elif platform.system() == "Darwin":
    lpm_lib = np.ctypeslib.load_library("liblpm_lib.dylib", "bin/")

np.random.seed(111)

curr_dir = os.getcwd()
model_lens = range(3, 10)
n_genes = 25
n_patients = 1000
fbp = 0.05
bgp = 0.05
n_reps = 100

_n_genes = ctypes.c_uint(n_genes)
_n_patients = ctypes.c_uint(n_patients)
_fbp = ctypes.c_double(fbp)
_bgp = ctypes.c_double(bgp)

for model_len in model_lens:
    _model_len = ctypes.c_uint(model_len)
    output_path = curr_dir + "/data/model_selection/model" + str(model_len)
    for rep in range(n_reps):
        dest = output_path + "/rep" + str(rep)
        if not os.path.exists(dest):
            os.makedirs(dest)
        _dest = ctypes.create_string_buffer(dest.encode())
        seed = np.random.randint(low=0, high=1000000)
        _seed = ctypes.c_long(seed)
        lpm_lib.generate_data(_seed, _dest, _model_len, _n_genes, _n_patients, _fbp, _bgp)

