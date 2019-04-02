import numpy as np 
from scipy.special import logsumexp
import ctypes
import os

np.random.seed(111)

lpm_lib = np.ctypeslib.load_library("liblpm_lib.dylib", "bin/")

lpm_lib.compute_likelihood.restype = ctypes.c_double

lpm_lib.run_smc.restype = ctypes.c_double
lpm_lib.run_smc.argtypes = [ctypes.c_long, ctypes.c_char_p, ctypes.c_uint, ctypes.c_uint,
                            ctypes.c_uint, ctypes.c_uint, ctypes.c_bool, ctypes.c_double,
                            ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_uint),
                            ctypes.POINTER(ctypes.c_double)]



curr_dir = os.getcwd()

seed = 1
n_mc_samples = 50
n_particles = 200
n_genes = 25
n_smc_iter = 2*n_genes
n_kernel_iter = 1
has_passenger = False
swap_prob = 0.2
error_prob_max = 0.1
n_threads = 8
true_model_len = 5
rep = 26

_seed = ctypes.c_long(seed)
_n_mc_samples = ctypes.c_uint(n_mc_samples)
_n_particles = ctypes.c_uint(n_particles)
_n_smc_iter = ctypes.c_uint(n_smc_iter)
_n_kernel_iter = ctypes.c_uint(n_kernel_iter)
_has_passenger = ctypes.c_bool(has_passenger)
_swap_prob = ctypes.c_double(swap_prob)
_n_threads = ctypes.c_uint(n_threads)
_n_genes = ctypes.c_uint(n_genes)
_true_model_len = ctypes.c_uint(true_model_len)

data_path = curr_dir + "/data/raphael/25genes/error0.05/rep" + str(rep) + "/"
input_path = data_path + "matrix.csv"
true_pathway_file = data_path + "generative_mem_mat.csv"
true_param_file = data_path + "parameters.csv"

_input_path = ctypes.create_string_buffer(input_path.encode())

# bgps and fbps
bgps = np.random.uniform(low=0.0, high=error_prob_max, size=n_mc_samples)
fbps = np.copy(bgps)

_bgps = np.ctypeslib.as_ctypes(bgps)
_fbps = np.ctypeslib.as_ctypes(fbps)

# compute the likelihood under the true state at the sampled parameters
pathway_matrix = np.genfromtxt(true_pathway_file, delimiter=',')
true_pathway = np.int32(np.argmax(pathway_matrix, axis=1))
_true_pathway = np.ctypeslib.as_ctypes(true_pathway)

true_param_file = data_path + "parameters.csv"
true_params = np.genfromtxt(true_param_file, delimiter=',')

_particles = (ctypes.c_uint * (n_particles * n_genes))()
_log_weights = (ctypes.c_double * n_particles)()
_fbp = ctypes.c_double(true_params[0])
_bgp = ctypes.c_double(true_params[1])
logZ = lpm_lib.run_smc(_seed, _input_path, _true_model_len, _n_particles, _n_smc_iter, 
                        _n_kernel_iter, _has_passenger, _swap_prob, _fbp, _bgp, 
                        _particles, _log_weights)

# check if the true state has been found in the samples
particles = np.reshape(np.ctypeslib.as_array(_particles), newshape=[n_particles, n_genes])
log_weights = np.reshape(np.ctypeslib.as_array(_log_weights), newshape=[n_particles, 1])
log_norm = logsumexp(log_weights)
for i in range(n_particles):
    particle = particles[i,:]
    _particle = np.ctypeslib.as_ctypes(np.int32(particle))
    a = np.sum(np.abs(particle - true_pathway))
    b = np.exp(log_weights[i] - log_norm)
    print(particles[i,:])
    print(a)
    print(b)
    log_lik = lpm_lib.compute_likelihood(_input_path, _particle, _true_model_len, 
                                            _n_genes, _has_passenger, _fbp, _bgp)
    print("loglik: " + str(log_lik))

print(true_pathway)
log_lik = lpm_lib.compute_likelihood(_input_path, _true_pathway, _true_model_len, 
                                            _n_genes, _has_passenger, _fbp, _bgp)
print("loglik@truth: " + str(log_lik))
