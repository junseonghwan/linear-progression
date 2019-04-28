import numpy as np 
from scipy.special import logsumexp
import ctypes
import os
import platform

if platform.system() == "Linux":
    lpm_lib = np.ctypeslib.load_library("liblpm_lib.so", "bin/")    
elif platform.system() == "Darwin":
    lpm_lib = np.ctypeslib.load_library("liblpm_lib.dylib", "bin/")

lpm_lib.compute_likelihood.restype = ctypes.c_double

np.random.seed(111)

curr_dir = os.getcwd()

seed = 1
n_genes = 100
n_pg_iter = 50
n_particles = 200
n_smc_iter = n_genes
n_kernel_iter = 5
n_mh_w_gibbs_iter = 10
has_passenger = True
swap_prob = 0.1
fbp_max = 0
bgp_max = 0.2
n_threads = 8
true_model_len = 5
n_patients = 1000
mh_proposal_sd = 0.05
rep = 0

_seed = ctypes.c_long(seed)
_n_genes = ctypes.c_uint(n_genes)
_n_pg_iter = ctypes.c_uint(n_pg_iter)
_n_particles = ctypes.c_uint(n_particles)
_n_smc_iter = ctypes.c_uint(n_smc_iter)
_n_kernel_iter = ctypes.c_uint(n_kernel_iter)
_n_mh_w_gibbs_iter = ctypes.c_uint(n_mh_w_gibbs_iter)
_has_passenger = ctypes.c_bool(has_passenger)
_swap_prob = ctypes.c_double(swap_prob)
_fbp_max = ctypes.c_double(fbp_max)
_bgp_max = ctypes.c_double(bgp_max)
_n_threads = ctypes.c_uint(n_threads)
_mh_proposal_sd = ctypes.c_double(mh_proposal_sd)
_true_model_len = ctypes.c_uint(true_model_len)

data_path = curr_dir + "/data/Experiment1/With_passengers/Uniform/error0.001/rep" + str(rep) + "/"
input_path = data_path + "matrix.csv"
output_path = data_path + "pg/"
true_pathway_file = data_path + "generative_mem_mat.csv"
true_param_file = data_path + "parameters.csv"

_input_path = ctypes.create_string_buffer(input_path.encode())
_output_path = ctypes.create_string_buffer(output_path.encode())

if not os.path.exists(output_path):
    os.makedirs(output_path)

lpm_lib.run_pg(_seed, _input_path, _output_path, _true_model_len, _n_pg_iter, 
                _n_particles, _n_smc_iter, _n_kernel_iter, _n_mh_w_gibbs_iter,
                _has_passenger, _swap_prob, _fbp_max, _bgp_max, _mh_proposal_sd)

# n_smc_iter = ctypes.c_uint(100)
# fbp = ctypes.c_double(0.01)
# bgp = ctypes.c_double(0.01)
# particles = (ctypes.c_uint * (100 * 25))()
# log_weights = (ctypes.c_double * 100)()
