
import numpy as np

# Determine the dataset you want to test:

error = 0.1
rep = 10
scenario = 'Decreasing'
passengers_status = 'Without_passengers'

dir = './Experiment1/'+passengers_status+'/'+scenario+'/error{}/rep{}/'.format(error, rep)

# Testing error rates:
a = np.genfromtxt(dir + 'clean_matrix.csv', delimiter=',')
b = np.genfromtxt(dir +'matrix.csv', delimiter=',')

number_of_errors_in_patients = np.sum(a != b, axis=1)
number_of_errors_in_genes = np.sum(a != b, axis=0)
total_number_of_errors = np.sum(number_of_errors_in_patients)
""" print('number_of_errors_in_patients')
print(number_of_errors_in_patients)
print('number_of_errors_in_genes')
print(number_of_errors_in_genes)
print('total_number_of_errors')
print(total_number_of_errors) """
if passengers_status == 'Without_passengers':
    error_rate = total_number_of_errors/(500*30)
else:
    error_rate = total_number_of_errors/(500*100)
print('Error_rate is {}'.format(error_rate))

# Testing pw capacity:
mem_mat = np.genfromtxt(dir + 'generative_mem_mat.csv', delimiter=',')
pathways_capacity = np.sum(mem_mat, axis=0)
print(pathways_capacity)

# Testing stages dist:
stages = np.genfromtxt(dir + 'patients_stages.csv', delimiter=',')
(unique_stages, counts)=np.unique(stages, return_counts=True)
print('observed stages:', unique_stages)
print('number of patients in each stage:', counts)