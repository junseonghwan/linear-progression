
import numpy as np

# Determine the dataset you want to test:

rep = 10
n_patients = 200
passengers_status = 'Without_passengers'

dir = './Experiment3/'+passengers_status+'/'+str(n_patients)+'/error0.1/rep{}/'.format(rep)

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
    error_rate = total_number_of_errors/(n_patients*30)
else:
    error_rate = total_number_of_errors/(n_patients*100)
print('Error_rate is {}'.format(error_rate))

print('Number of patients is {}'.format(a.shape[0]))

# Testing pw capacity:
mem_mat = np.genfromtxt(dir + 'generative_mem_mat.csv', delimiter=',')
pathways_capacity = np.sum(mem_mat, axis=0)
print(pathways_capacity)

# Testing stages dist:
stages = np.genfromtxt(dir + 'patients_stages.csv', delimiter=',')
(unique_stages, counts)=np.unique(stages, return_counts=True)
print('observed stages:', unique_stages)
print('number of patients in each stage:', counts)