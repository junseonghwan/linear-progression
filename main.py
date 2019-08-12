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
configs = functions.parse_config(config_file)

seed = int(configs["seed"])
model_lens = configs["model_len"].split("-")
if len(model_lens) == 2:
    model_len_begin = int(model_lens[0])
    model_len_end = int(model_lens[1])
else:
    model_len_begin = int(model_lens[0])
    model_len_end = model_len_begin
n_mcmc_iter = int(configs["n_mcmc_iter"])
n_mh_w_gibbs_iter = int(configs["n_mh_w_gibbs_iter"])
has_passenger = bool(configs["has_passenger"].lower() == "true")
swap_prob = float(configs["swap_prob"])
fbp_max = float(configs["fbp_max"])
bgp_max = float(configs["bgp_max"])
mh_proposal_sd = float(configs["mh_proposal_sd"])
data_path = os.path.abspath(configs["data_path"])
output_path = os.path.abspath(configs["output_path"])
thinning_interval = int(configs["thinning_interval"])
burn_in = int(configs["burn_in"])
prior_passenger_prob = float(configs["prior_passenger_prob"])

data_file = data_path
curr_time = str(time.time()).encode("utf-8")
dir_name = hashlib.sha1(curr_time).hexdigest()
output_path = output_path + "/" + dir_name
if not os.path.exists(output_path):
    os.makedirs(output_path)
for model_len in range(model_len_begin, model_len_end+1):
    output_file = output_path + "/model" + str(model_len)
    if not os.path.exists(output_file):
        os.makedirs(output_file)
    cmd = ["./bin/lpm"]
    functions.add_cmd(cmd, "-d", data_file)
    functions.add_cmd(cmd, "-o", output_file)
    functions.add_cmd(cmd, "-s", seed)
    functions.add_cmd(cmd, "-l", model_len)
    functions.add_cmd(cmd, "-p", n_mcmc_iter)
    functions.add_cmd(cmd, "-G", n_mh_w_gibbs_iter)
    functions.add_cmd(cmd, "-f", fbp_max)
    functions.add_cmd(cmd, "-b", bgp_max)
    functions.add_cmd(cmd, "--has_passenger", has_passenger)
    functions.add_cmd(cmd, "--swap_prob", swap_prob)
    functions.add_cmd(cmd, "--mh_proposal_sd", mh_proposal_sd)
    functions.add_cmd(cmd, "--burn_in", burn_in)
    functions.add_cmd(cmd, "--thinning_interval", thinning_interval)
    functions.add_cmd(cmd, "--prior_passenger_prob", prior_passenger_prob)

    # run LPM
    start = time.time()
    subprocess.run(cmd)
    # copy config.txt
    copyfile(config_file, output_path + "/" + config_file)
    end = time.time()
    print(str(end - start) + " seconds.")
