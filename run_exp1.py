from scipy.special import logsumexp
import numpy as np 
import ctypes
import os
import platform
import sys
import functions

config_file = sys.argv[1]
rep_begin = int(sys.argv[2])
rep_end = int(sys.argv[3])

key_list = ["seed", "model_len", "n_mc_samples", "n_pg_iter", "n_particles", "n_smc_iter", "n_kernel_iter", 
            "n_mh_w_gibbs_iter", "has_passenger", "swap_prob", "fbp_max", "bgp_max", "mh_proposal_sd", 
            "n_threads", "data_path", "use_lik_tempering"]

configs = functions.parse_config(config_file)
if not functions.check_configs(configs, key_list):
    sys.exit(-1)

if platform.system() == "Linux":
    lpm_lib = np.ctypeslib.load_library("liblpm_lib.so", "bin/")
elif platform.system() == "Darwin":
    lpm_lib = np.ctypeslib.load_library("liblpm_lib.dylib", "bin/")
lpm_lib.model_selection.restype = ctypes.c_double

seed = int(configs["seed"])
model_lens = configs["model_len"].split("-")
if len(model_lens) == 2:
    model_len_begin = int(model_lens[0])
    model_len_end = int(model_lens[1])
else:
    model_len_begin = 2
    model_len_end = 10
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
use_lik_tempering = bool(configs["use_lik_tempering"])
true_model_len = 5

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
_use_lik_tempering = ctypes.c_bool(use_lik_tempering)
_true_model_len = ctypes.c_uint(true_model_len)

# configs["data_path"] contains directories labelled: rep[0-9]+
for rep in range(rep_begin, rep_end+1):
    rep_path = data_path + "/rep" + str(rep) + "/"
    input_path = rep_path + "matrix.csv"
    model_selection_output_path = rep_path + "model_selection/"
    pg_output_path = rep_path + "pg_true_len/"
    auto_pg_output_path = rep_path + "pg_auto/"
    _input_path = ctypes.create_string_buffer(input_path.encode())

    # 1. run model selection
    bgps = functions.generate_stratified_samples(bgp_max, n_mc_samples)
    if fbp_max > 0.0:
        fbps = functions.generate_stratified_samples(fbp_max, n_mc_samples)
    else:
        fbps = np.copy(bgps)

    _bgps = np.ctypeslib.as_ctypes(bgps)
    _fbps = np.ctypeslib.as_ctypes(fbps)

    _log_marginals = (ctypes.c_double * n_mc_samples)()
    _log_marginals_smc = (ctypes.c_double * n_mc_samples)()

    fhats = np.zeros(shape=(model_len_end - model_len_begin + 1,2))
    log_marginals_matrix = []
    for model_len in range(model_len_begin, model_len_end + 1):
        _model_len = ctypes.c_uint(model_len)
        fhat = lpm_lib.model_selection(seed, _input_path, _model_len, _n_mc_samples, 
                                        _n_particles, _n_smc_iter, _n_kernel_iter, 
                                        _has_passenger, _swap_prob, _fbps, _bgps,
                                        _n_threads, _log_marginals, _log_marginals_smc, 
                                        _use_lik_tempering)
        fhats[model_len-model_len_begin][0] = model_len
        fhats[model_len-model_len_begin][1] = fhat
        log_marginal = np.ctypeslib.as_array(_log_marginals)
        log_marginals_smc = np.ctypeslib.as_array(_log_marginals_smc)
        log_marginals_matrix.append(np.column_stack((np.repeat(model_len, n_mc_samples), fbps, bgps, log_marginal, log_marginals_smc)))

    if not os.path.exists(model_selection_output_path):
        os.makedirs(model_selection_output_path)

    fhat_file = model_selection_output_path + "fhat.csv"
    np.savetxt(fname=fhat_file, X=fhats, fmt="%d,%f", header="Model,f")

    dat = np.concatenate(log_marginals_matrix, axis=0)
    log_marginals_file = model_selection_output_path + "log_marginals.csv"
    np.savetxt(fname=log_marginals_file, X=dat, fmt="%d,%f,%f,%f, %f", header="Model,FBP,BGP,MarginalLogLikSum,MarginalLogLikSMC")

    # 2. run PG using the best model len
    best_model_len = int(fhats[np.argmax(fhats[:,1]),0])
    _best_model_len = ctypes.c_uint(best_model_len)
    _output_path = ctypes.create_string_buffer(auto_pg_output_path.encode())

    if not os.path.exists(auto_pg_output_path):
        os.makedirs(auto_pg_output_path)

    lpm_lib.run_pg(_seed, _input_path, _output_path, _best_model_len, _n_pg_iter, 
                    _n_particles, _n_smc_iter, _n_kernel_iter, _n_mh_w_gibbs_iter,
                    _has_passenger, _swap_prob, _fbp_max, _bgp_max, _mh_proposal_sd)

    # 3. run PG using the true model len if necessary
    if best_model_len != true_model_len:
        _output_path = ctypes.create_string_buffer(pg_output_path.encode())
        if not os.path.exists(pg_output_path):
            os.makedirs(pg_output_path)

        lpm_lib.run_pg(_seed, _input_path, _output_path, _true_model_len, _n_pg_iter, 
                        _n_particles, _n_smc_iter, _n_kernel_iter, _n_mh_w_gibbs_iter,
                        _has_passenger, _swap_prob, _fbp_max, _bgp_max, _mh_proposal_sd)
