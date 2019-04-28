# Python script to read in the config file and run LPM
import numpy as np 
from scipy.special import logsumexp
import ctypes
import os
import platform
import sys
import functions

if platform.system() == "Linux":
    lpm_lib = np.ctypeslib.load_library("liblpm_lib.so", "bin/")
elif platform.system() == "Darwin":
    lpm_lib = np.ctypeslib.load_library("liblpm_lib.dylib", "bin/")

lpm_lib.model_selection.restype = ctypes.c_double

key_list = ["seed", "model_len", "n_mc_samples", "n_pg_iter", "n_particles", "n_smc_iter", "n_kernel_iter", 
            "n_mh_w_gibbs_iter", "has_passenger", "swap_prob", "fbp_max", "bgp_max", "mh_proposal_sd", 
            "n_threads", "data_path", "output_path", "use_lik_tempering"]

config_file = sys.argv[1]
configs = functions.parse_config(config_file)
if not functions.check_configs(configs, key_list):
    sys.exit(-1)

seed = int(configs["seed"])
model_lens = configs["model_len"].split("-")
if len(model_lens) == 2:
    model_len_begin = int(model_lens[0])
    model_len_end = int(model_lens[1])
else:
    model_len_begin = int(model_lens[0])
    model_len_end = model_len_begin
n_mc_samples = int(configs["n_mc_samples"])
n_pg_iter = int(configs["n_pg_iter"])
n_particles = int(configs["n_particles"])
n_smc_iter = int(configs["n_smc_iter"])
n_kernel_iter = int(configs["n_kernel_iter"])
n_mh_w_gibbs_iter = int(configs["n_mh_w_gibbs_iter"])
has_passenger = bool(configs["has_passenger"])
swap_prob = float(configs["swap_prob"])
fbp_max = float(configs["fbp_max"])
bgp_max = float(configs["bgp_max"])
mh_proposal_sd = float(configs["mh_proposal_sd"])
n_threads = int(configs["n_threads"])
data_path = os.path.abspath(configs["data_path"])
output_path = os.path.abspath(configs["output_path"])
use_lik_tempering = bool(configs["use_lik_tempering"])

_seed = ctypes.c_long(seed)
_n_mc_samples = ctypes.c_uint(n_mc_samples)
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
_data_path = ctypes.create_string_buffer(data_path.encode())
_use_lik_tempering = ctypes.c_bool(use_lik_tempering)

# run model selection to choose the model_len then, run PG
if model_len_begin < model_len_end:
    # stratification
    bgps = functions.generate_stratified_samples(bgp_max, n_mc_samples)
    if fbp_max > 0.0:
        fbps = functions.generate_stratified_samples(fbp_max, n_mc_samples)
    else:
        fbps = np.copy(bgps)

    print(bgps)
    _bgps = np.ctypeslib.as_ctypes(bgps)
    _fbps = np.ctypeslib.as_ctypes(fbps)

    _log_marginals = (ctypes.c_double * n_mc_samples)()
    _log_marginals_smc = (ctypes.c_double * n_mc_samples)()

    fhats = np.zeros(shape=(model_len_end - model_len_begin + 1,2))
    log_marginals_matrix = []
    for model_len in range(model_len_begin, model_len_end + 1):
        _model_len = ctypes.c_uint(model_len)
        fhat = lpm_lib.model_selection(seed, _data_path, _model_len, _n_mc_samples, 
                                        _n_particles, _n_smc_iter, _n_kernel_iter, 
                                        _has_passenger, _swap_prob, _fbps, _bgps,
                                        _n_threads, _log_marginals, _log_marginals_smc, 
                                        _use_lik_tempering)
        fhats[model_len-model_len_begin][0] = model_len
        fhats[model_len-model_len_begin][1] = fhat
        log_marginal = np.ctypeslib.as_array(_log_marginals)
        log_marginals_smc = np.ctypeslib.as_array(_log_marginals_smc)
        log_marginals_matrix.append(np.column_stack((np.repeat(model_len, n_mc_samples), fbps, bgps, log_marginal, log_marginals_smc)))

    output_path_model_selection = output_path + "/model_selection/"
    if not os.path.exists(output_path_model_selection):
        os.makedirs(output_path_model_selection)

    fhat_file = output_path_model_selection + "fhat.csv"
    np.savetxt(fname=fhat_file, X=fhats, fmt="%d,%f", header="Model,f")

    dat = np.concatenate(log_marginals_matrix, axis=0)
    log_marginals_file = output_path_model_selection + "log_marginals.csv"
    np.savetxt(fname=log_marginals_file, X=dat, fmt="%d,%f,%f,%f, %f", header="Model,FBP,BGP,MarginalLogLikSum,MarginalLogLikSMC")

    model_len = int(fhats[np.argmax(fhats[:,1]),0])
else:
    model_len = model_len_begin

print(model_len)
_model_len = ctypes.c_uint(model_len)

# run PG
output_path_pg = output_path + "/pg/"
_output_path_pg = ctypes.create_string_buffer(output_path_pg.encode())
if not os.path.exists(output_path_pg):
    os.makedirs(output_path_pg)

lpm_lib.run_pg(_seed, _data_path, _output_path_pg, _model_len, _n_pg_iter, 
                _n_particles, _n_smc_iter, _n_kernel_iter, _n_mh_w_gibbs_iter,
                _has_passenger, _swap_prob, _fbp_max, _bgp_max, _mh_proposal_sd)
