import numpy as np 
import ctypes
import os

lpm_lib = np.ctypeslib.load_library("liblpm_lib.dylib", "build/src/")

curr_dir = os.getcwd()
data_path = curr_dir + "/data/raphael/25genes/error0.01/rep1/"
input_path = data_path + "matrix.csv"
output_path = data_path + "pg/"

# check output_path exists
if not os.path.exists(output_path):
    os.makedirs(output_path)

seed = ctypes.c_long(1)
_input_path = ctypes.create_string_buffer(input_path.encode())
_output_path = ctypes.create_string_buffer(output_path.encode())
model_len = ctypes.c_uint(5)
n_pg_iter = ctypes.c_uint(10)
n_particles = ctypes.c_uint(100)
n_smc_iter = ctypes.c_uint(25)
n_kernel_iter = ctypes.c_uint(1)
n_mh_w_gibbs_iter = ctypes.c_uint(20)
has_passenger = ctypes.c_bool(False)
swap_prob = ctypes.c_double(0.2)
fbp_max = ctypes.c_double(0.0)
bgp_max = ctypes.c_double(0.3)

# uncomment below to run PG
# lpm_lib.run_pg(seed, _input_path, _output_path, model_len, n_pg_iter, 
#                 n_particles, n_smc_iter, n_kernel_iter, n_mh_w_gibbs_iter,
#                 has_passenger, swap_prob, fbp_max, bgp_max)

n_smc_iter = ctypes.c_uint(100)
fbp = ctypes.c_double(0.01)
bgp = ctypes.c_double(0.01)
particles = (ctypes.c_uint * (100 * 25))()
log_weights = (ctypes.c_double * 100)()
lpm_lib.run_smc.restype = ctypes.c_double
# SMC returns logZ, as well as particles and log weights (the latter two are passed into the function)
lpm_lib.run_smc.argtypes = [ctypes.c_long, ctypes.c_char_p, ctypes.c_uint, ctypes.c_uint,
                            ctypes.c_uint, ctypes.c_uint, ctypes.c_bool, ctypes.c_double,
                            ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_uint),
                            ctypes.POINTER(ctypes.c_double)]
logZ = lpm_lib.run_smc(seed, _input_path, model_len, n_particles, n_smc_iter, n_kernel_iter, 
                has_passenger, swap_prob, fbp, bgp, particles, log_weights)
print(logZ)

particles_arr = np.reshape(np.ctypeslib.as_array(particles), newshape=[100, 25])
print(particles_arr.shape)
for i in range(100):
    print("%f" % log_weights[i])
    for j in range(25):
        print("%d " % particles_arr[i,j], end='')
    print("")
