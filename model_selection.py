import numpy as np 
import ctypes
import os

np.random.seed(123)

lpm_lib = np.ctypeslib.load_library("liblpm_lib.dylib", "bin/")
curr_dir = os.getcwd()
data_path = curr_dir + "/data/raphael/25genes/error0.05/rep1/"
output_path = data_path + "model_selection/"
input_path = data_path + "matrix.csv"
true_pathway_file = data_path + "generative_mem_mat.csv"
true_param_file = data_path + "parameters.csv"

if not os.path.exists(output_path):
    os.makedirs(output_path)

seed = 1
n_mc_samples = 50
n_particles = 100
n_genes = 25
n_smc_iter = n_genes
n_kernel_iter = 1
has_passenger = False
swap_prob = 0.2
error_prob_max = 0.1
n_threads = 8
model_begin = 2
model_end = 10

_seed = ctypes.c_long(seed)
_input_path = ctypes.create_string_buffer(input_path.encode())
_n_mc_samples = ctypes.c_uint(n_mc_samples)
_n_particles = ctypes.c_uint(n_particles)
_n_smc_iter = ctypes.c_uint(n_smc_iter)
_n_kernel_iter = ctypes.c_uint(n_kernel_iter)
_has_passenger = ctypes.c_bool(has_passenger)
_swap_prob = ctypes.c_double(swap_prob)
_n_threads = ctypes.c_uint(n_threads)
_n_genes = ctypes.c_uint(n_genes)
_log_marginals = (ctypes.c_double * n_mc_samples)()

# compute the likelihood under the true state at the sampled parameters
pathway_matrix = np.genfromtxt(true_pathway_file, delimiter=',')
true_pathway = np.int32(np.argmax(pathway_matrix, axis=1))
_true_pathway = np.ctypeslib.as_ctypes(true_pathway)

true_model_len = pathway_matrix.shape[1]
_true_model_len = ctypes.c_uint(true_model_len)

true_params = np.genfromtxt(true_param_file, delimiter=',')
_bgp = ctypes.c_double(true_params[0])
_fbp = ctypes.c_double(true_params[1])

# generate fbps and bgps from the prior distribution
bgps = np.random.uniform(low=0.0, high=error_prob_max, size=n_mc_samples)
fbps = np.copy(bgps)

_bgps = np.ctypeslib.as_ctypes(bgps)
_fbps = np.ctypeslib.as_ctypes(fbps)

# compute the log likelihood under the sampled parameters
lpm_lib.compute_likelihood.restype = ctypes.c_double
true_state_log_liks = np.zeros(shape=(n_mc_samples,3))
for i in range(n_mc_samples):
    _bgp = ctypes.c_double(bgps[i])
    _fbp = ctypes.c_double(fbps[i])
    log_lik = lpm_lib.compute_likelihood(_input_path, _true_pathway, _true_model_len, 
                                        _n_genes, _has_passenger, _fbp, _bgp)
    true_state_log_liks[i][0] = fbps[i]
    true_state_log_liks[i][1] = bgps[i]
    true_state_log_liks[i][2] = log_lik

# output the log_liks under true state
log_lik_at_truth_file = output_path + "log_lik_at_truth.csv"
np.savetxt(fname=log_lik_at_truth_file, X=true_state_log_liks, fmt="%f", delimiter=',',header="FBP,BGP,LogLik")

# run model selection experiments output \hat{f} and 
lpm_lib.model_selection.restype = ctypes.c_double
fhats = np.zeros(shape=(model_end - model_begin + 1,2))
log_marginals_matrix = []
for model_len in range(model_begin, model_end + 1):
    _model_len = ctypes.c_uint(model_len)
    fhat = lpm_lib.model_selection(seed, _input_path, _model_len, _n_mc_samples, 
                                    _n_particles, _n_smc_iter, _n_kernel_iter, 
                                    _has_passenger, _swap_prob, _fbps, _bgps,
                                    _n_threads, _log_marginals)
    print("fhat: " + str(fhat))
    fhats[model_len-model_begin][0] = model_len
    fhats[model_len-model_begin][1] = fhat
    log_marginal = np.ctypeslib.as_array(_log_marginals)
    log_marginals_matrix.append(np.column_stack((np.repeat(model_len, n_mc_samples), fbps, bgps, log_marginal)))

fhat_file = output_path + "fhat.csv"
np.savetxt(fname=fhat_file, X=fhats, fmt="%d,%f", header="Model,f")

dat = np.concatenate(log_marginals_matrix, axis=0)
log_marginals_file = output_path + "log_marginals.csv"
np.savetxt(fname=log_marginals_file, X=dat, fmt="%d,%f,%f,%f", header="Model,FBP,BGP,MarginalLogLik")