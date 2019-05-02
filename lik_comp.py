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

curr_dir = os.getcwd()

rep = 2
data_path = curr_dir + "/data/Experiment1/With_passengers/Uniform/error0.001/rep" + str(rep) + "/"
input_path = data_path + "matrix.csv"
_input_path = ctypes.create_string_buffer(input_path.encode())

true_pathway_file = data_path + "generative_mem_mat.csv"
pathway_matrix = np.genfromtxt(true_pathway_file, delimiter=',')
true_pathway = np.int32(np.argmax(pathway_matrix, 1))
_true_pathway = np.ctypeslib.as_ctypes(true_pathway)

model_len = 5
_model_len = ctypes.c_uint(model_len)

n_genes = 100
_n_genes = ctypes.c_uint(n_genes)

has_passenger = True
_has_passenger = ctypes.c_bool(has_passenger)

bgp = 0.001
fbp = 0.001
_bgp = ctypes.c_double(bgp)
_fbp = ctypes.c_double(fbp)

log_lik = lpm_lib.compute_likelihood(_input_path, _true_pathway, _model_len, _n_genes, _has_passenger, _fbp, _bgp)
print(log_lik)
