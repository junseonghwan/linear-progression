
import os
import numpy as np

from dataset import Dataset


def synthetic_data_construction_fer(number_of_datasets=100,
                                    PE=True,
                                    n_genes=100,
                                    model_length=6,
                                    pw_capacity=[6, 6, 6, 6, 6],
                                    error_param_list=[0.001, 0.01, 0.1],
                                    n_patients=500,
                                    number_of_passenger_genes=None,
                                    output_path='./Datasets/'):

    for error_param in error_param_list:
        if 'error{}'.format(error_param) not in os.listdir(output_path):
            os.mkdir(output_path+'error{}/'.format(error_param))
        for i in np.arange(number_of_datasets):
            if 'rep{}'.format(i) not in os.listdir(output_path+'error{}/'.format(error_param)):
                os.mkdir(output_path+'error{}/rep{}/'.format(error_param, i))
            dataset = Dataset.from_parameters(
                n_genes=n_genes,
                model_length=model_length,
                n_genes_in_pws=pw_capacity,
                generative_fbp=error_param,
                generative_bgp=error_param,
                n_patients=n_patients,
                PE=PE,
                number_of_passenger_genes=number_of_passenger_genes
                )
            dataset.save_data(dir=output_path+'error{}/rep{}/'.format(error_param, i))
    return()


def synthetic_data_construction_fen():
    # fixed error number
    number_of_datasets = 10

    if 'Datasets_fen' not in os.listdir():
        os.mkdir('./Datasets_fen/')

    PE = True
    n_genes = 100
    pw_capacity = [5, 5, 5, 5, 5]
    model_length = 6
    n_patients_list = [50, 100, 200, 500, 1000]
    error_per_patient_list = n_genes*np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.2])

    if '{}genes'.format(n_genes) not in os.listdir('./Datasets_fen'):
        os.mkdir('./Datasets_fen/{}genes'.format(n_genes))
    for n_patients in n_patients_list:
        if 'n_patients{}'.format(n_patients) not in os.listdir('./Datasets_fen/{}genes'.format(n_genes)):
            os.mkdir('./Datasets_fen/{}genes/n_patients{}'.format(n_genes, n_patients))
        n_errors_list = (n_patients*error_per_patient_list).astype(int)
        for n_errors in n_errors_list:
            if 'n_errors{}'.format(n_errors) not in os.listdir('./Datasets_fen/{}genes/n_patients{}'.format(n_genes, n_patients)):
                    os.mkdir('./Datasets_fen/{}genes/n_patients{}/n_errors{}'.format(n_genes, n_patients, n_errors))
            for i in np.arange(number_of_datasets):
                if 'rep{}'.format(i) not in os.listdir('./Datasets_fen/{}genes/n_patients{}/n_errors{}'.format(n_genes, n_patients, n_errors)):
                    os.mkdir('./Datasets_fen/{}genes/n_patients{}/n_errors{}/rep{}'.format(n_genes, n_patients, n_errors, i))
                dataset = Dataset.from_parameters(
                    n_genes=n_genes,
                    model_length=model_length,
                    n_genes_in_pws=pw_capacity,
                    number_of_errors=n_errors,
                    n_patients=n_patients,
                    PE=PE
                    )
                dataset.save_data(dir='./Datasets_fen/{}genes/n_patients{}/n_errors{}/rep{}/'.format(n_genes, n_patients, n_errors, i))
    return()


if __name__ == '__main__':

    if 'Experiment1' not in os.listdir():
        os.mkdir('./Experiment1/')

    ## Experiment1 without passengers:
    if 'Without_passengers' not in os.listdir('./Experiment1/'):
        os.mkdir('./Experiment1/Without_passengers/')
    # Scenario1: Uniform pathway capacity
    if 'Uniform' not in os.listdir('./Experiment1/Without_passengers/'):
        os.mkdir('./Experiment1/Without_passengers/Uniform/')
    synthetic_data_construction_fer(
        number_of_datasets=100,
        PE=False,
        n_genes=30,
        model_length=5,
        pw_capacity=[6, 6, 6, 6],
        error_param_list=[0.001, 0.01, 0.1],
        n_patients=500,
        output_path='./Experiment1/Without_passengers/Uniform/'
    )
    # Scenario2: Increasing pathway capacity
    if 'Increasing' not in os.listdir('./Experiment1/Without_passengers/'):
        os.mkdir('./Experiment1/Without_passengers/Increasing/')
    synthetic_data_construction_fer(
        number_of_datasets=100,
        PE=False,
        n_genes=30,
        model_length=5,
        pw_capacity=[2, 4, 6, 8],
        error_param_list=[0.001, 0.01, 0.1],
        n_patients=500,
        output_path='./Experiment1/Without_passengers/Increasing/'
    )
    # Scenario3: Decreasing pathway capacity
    if 'Decreasing' not in os.listdir('./Experiment1/Without_passengers/'):
        os.mkdir('./Experiment1/Without_passengers/Decreasing/')
    synthetic_data_construction_fer(
        number_of_datasets=100,
        PE=False,
        n_genes=30,
        model_length=5,
        pw_capacity=[10, 8, 6, 4],
        error_param_list=[0.001, 0.01, 0.1],
        n_patients=500,
        output_path='./Experiment1/Without_passengers/Decreasing/'
    )
    # Scenario4: Random split
    if 'Random' not in os.listdir('./Experiment1/Without_passengers/'):
        os.mkdir('./Experiment1/Without_passengers/Random/')
    synthetic_data_construction_fer(
        number_of_datasets=100,
        PE=False,
        n_genes=30,
        model_length=5,
        pw_capacity=None,
        error_param_list=[0.001, 0.01, 0.1],
        n_patients=500,
        output_path='./Experiment1/Without_passengers/Random/'
    )

    ## Experiment1 with passengers:
    if 'With_passengers' not in os.listdir('./Experiment1/'):
        os.mkdir('./Experiment1/With_passengers/')
    # Scenario1: Uniform pathway capacity
    if 'Uniform' not in os.listdir('./Experiment1/With_passengers/'):
        os.mkdir('./Experiment1/With_passengers/Uniform/')
    synthetic_data_construction_fer(
        number_of_datasets=100,
        PE=True,
        n_genes=100,
        model_length=6,
        pw_capacity=[6, 6, 6, 6, 6],
        error_param_list=[0.001, 0.01, 0.1],
        n_patients=500,
        output_path='./Experiment1/With_passengers/Uniform/'
    )
    # Scenario2: Increasing pathway capacity
    if 'Increasing' not in os.listdir('./Experiment1/With_passengers/'):
        os.mkdir('./Experiment1/With_passengers/Increasing/')
    synthetic_data_construction_fer(
        number_of_datasets=100,
        PE=True,
        n_genes=100,
        model_length=6,
        pw_capacity=[2, 4, 6, 8, 10],
        error_param_list=[0.001, 0.01, 0.1],
        n_patients=500,
        output_path='./Experiment1/With_passengers/Increasing/'
    )
    # Scenario3: Decreasing pathway capacity
    if 'Decreasing' not in os.listdir('./Experiment1/With_passengers/'):
        os.mkdir('./Experiment1/With_passengers/Decreasing/')
    synthetic_data_construction_fer(
        number_of_datasets=100,
        PE=True,
        n_genes=100,
        model_length=6,
        pw_capacity=[10, 8, 6, 4, 2],
        error_param_list=[0.001, 0.01, 0.1],
        n_patients=500,
        output_path='./Experiment1/With_passengers/Decreasing/'
    )
    # Scenario4: Random split
    if 'Random' not in os.listdir('./Experiment1/With_passengers/'):
        os.mkdir('./Experiment1/With_passengers/Random/')
    synthetic_data_construction_fer(
        number_of_datasets=100,
        PE=True,
        n_genes=100,
        model_length=6,
        number_of_passenger_genes=70,
        pw_capacity=None,
        error_param_list=[0.001, 0.01, 0.1],
        n_patients=500,
        output_path='./Experiment1/With_passengers/Random/'
    )

    print('done')
