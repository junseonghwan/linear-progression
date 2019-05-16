# Python script to read in the config file and run LPM
import hashlib
import glob
import os
import sys
import subprocess
import time
from shutil import copyfile
import functions
from numba import jit

config_file = sys.argv[1]
mode = sys.argv[2] # {model, mcmc}
rep_begin = int(sys.argv[3])
rep_end = int(sys.argv[4])
configs = functions.parse_config(config_file)

seed = int(configs["seed"])
model_lens = configs["model_len"].split("-")
if len(model_lens) == 2:
    model_len_begin = int(model_lens[0])
    model_len_end = int(model_lens[1])
else:
    model_len_begin = int(model_lens[0])
    model_len_end = model_len_begin
n_mc_samples = int(configs["n_mc_samples"])
n_mcmc_iter = int(configs["n_mcmc_iter"])
n_particles = int(configs["n_particles"])
n_smc_iter = int(configs["n_smc_iter"])
n_kernel_iter = int(configs["n_kernel_iter"])
n_mh_w_gibbs_iter = int(configs["n_mh_w_gibbs_iter"])
has_passenger = bool(configs["has_passenger"])
swap_prob = float(configs["swap_prob"])
fbp_max = float(configs["fbp_max"])
bgp_max = float(configs["bgp_max"])
mh_proposal_sd = float(configs["mh_proposal_sd"])
n_mc_jobs = int(configs["n_mc_jobs"])
n_smc_threads = int(configs["n_smc_threads"])
data_path = os.path.abspath(configs["data_path"])
use_lik_tempering = bool(configs["use_lik_tempering"])
thinning_interval = int(configs["thinning_interval"])

for rep in range(rep_begin, rep_end+1):
    data_file = data_path + "/rep" + str(rep) + "/matrix.csv"
    mode_output_path = data_path + "/rep" + str(rep) + "/" + mode
    curr_time = str(time.time()).encode("utf-8")
    dir_name = hashlib.sha1(curr_time).hexdigest()
    output_path = mode_output_path + "/" + dir_name
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for model_len in range(model_len_begin, model_len_end+1):
        if mode == "model":
            output_file = output_path + "/log_marginals" + str(model_len) + ".csv"
        elif mode == "mcmc":
            output_file = output_path + "/model" + str(model_len)
            if not os.path.exists(output_file):
                os.makedirs(output_file)
        cmd = ["./bin/lpm"]
        functions.add_cmd(cmd, "-d", data_file)
        functions.add_cmd(cmd, "-o", output_file)
        functions.add_cmd(cmd, "-m", mode)
        functions.add_cmd(cmd, "-s", seed)
        functions.add_cmd(cmd, "-l", model_len)
        functions.add_cmd(cmd, "-M", n_mc_samples)
        functions.add_cmd(cmd, "-p", n_mcmc_iter)
        functions.add_cmd(cmd, "-P", n_particles)
        functions.add_cmd(cmd, "-S", n_smc_iter)
        functions.add_cmd(cmd, "-k", n_kernel_iter)
        functions.add_cmd(cmd, "-G", n_mh_w_gibbs_iter)
        functions.add_cmd(cmd, "-t", n_mc_jobs)
        functions.add_cmd(cmd, "--n_smc_threads", n_smc_threads)
        functions.add_cmd(cmd, "-T", use_lik_tempering)
        functions.add_cmd(cmd, "-f", fbp_max)
        functions.add_cmd(cmd, "-b", bgp_max)
        functions.add_cmd(cmd, "--has_passenger", has_passenger)
        functions.add_cmd(cmd, "--swap_prob", swap_prob)
        functions.add_cmd(cmd, "--mh_proposal_sd", mh_proposal_sd)
        functions.add_cmd(cmd, "--thinning_interval", thinning_interval)

        # run LPM
        subprocess.run(cmd)
        # copy config.txt
        copyfile(config_file, output_path + "/" + config_file)
