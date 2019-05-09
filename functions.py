import numpy as np

def parse_config(config_file):
    f = open(config_file, "r")
    configs = dict()
    for line in f:
        row = line.split(":")
        configs[row[0].strip()] = row[1].strip()
    return configs

def check_configs(configs, key_list):
    for key in key_list:
        if key not in configs:
            print("Error: " + key + " not in the configuration file.")
            return False
    return True

def generate_stratified_samples(max_val, num_samples):
    ret_list = []
    strata = np.arange(0, max_val, max_val/num_samples)
    for i, val in enumerate(strata):
        begin = strata[i]
        if (i+1) == len(strata):
            end = max_val
        else:
            end = strata[i+1]
        ret_list.append(np.random.uniform(low=begin, high=end, size=1)[0])
    return np.asarray(ret_list)

def add_cmd(cmd, option, val):
    cmd.append(option)
    cmd.append(str(val))
