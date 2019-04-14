
import os
import copy
import numpy as np


class Dataset:

    def __init__(self, matrix, true_model_dict=None):

        self.matrix = matrix
        (self.n_patients, self.n_genes) = np.shape(matrix)

        if true_model_dict is not None:

            self.generative_mem_mat = true_model_dict['generative_mem_mat']
            self.generative_K = np.shape(self.generative_mem_mat)[1]
            self.patients_stages = true_model_dict['patients_stages']
            self.clean_matrix = true_model_dict['clean_matrix']

            if 'generative_bgp' in true_model_dict.keys():
                self.generative_bgp = true_model_dict['generative_bgp']
            if 'generative_fbp' in true_model_dict.keys():
                self.generative_fbp = true_model_dict['generative_fbp']
            if 'number_of_errors' in true_model_dict.keys():
                self.number_of_errors = true_model_dict['number_of_errors']

    @property
    def generative_pathways_model(self):
        n_genes = self.n_genes
        model_length = self.generative_K
        pathways_model = [[] for i in range(0, model_length)]
        for i in range(0, n_genes):
            for j in range(0, model_length):
                if self.generative_mem_mat[i, j] == 1:
                    pathways_model[j].append(i)
        return(pathways_model)

    def save_data(self, dir='./'):

        np.savetxt(dir + 'matrix.csv', self.matrix, fmt='%d', delimiter=',')

        if 'clean_matrix' in self.__dict__:
            np.savetxt(
                dir + 'clean_matrix.csv',
                self.clean_matrix,
                fmt='%d',
                delimiter=','
                )
        if 'generative_mem_mat' in self.__dict__:
            np.savetxt(
                dir + 'generative_mem_mat.csv',
                self.generative_mem_mat,
                fmt='%d',
                delimiter=','
                )
        if 'patients_stages' in self.__dict__:
            np.savetxt(
                dir + 'patients_stages.csv',
                self.patients_stages,
                fmt='%d',
                delimiter=','
                )
        if 'generative_bgp' in self.__dict__ and \
           'generative_fbp' in self.__dict__:
            parameters = np.array([self.generative_fbp, self.generative_bgp])
            np.savetxt(
                dir + 'parameters.csv',
                parameters,
                delimiter=','
                )
        if 'number_of_errors' in self.__dict__:
            np.savetxt(
                dir + 'parameters.csv',
                np.array([self.number_of_errors]),
                delimiter=','
                )
        return()

    @classmethod
    def from_csv(cls, path='./', n_patients=None):
        matrix = np.genfromtxt(path + '/matrix.csv', delimiter=',').astype(int)
        true_model_dict = {}
        if os.path.isfile(path + '/patients_stages.csv'):
            true_model_dict['patients_stages'] = np.genfromtxt(
                path + '/patients_stages.csv', delimiter=','
                ).astype(int)

        if os.path.isfile(path + '/generative_mem_mat.csv'):
            true_model_dict['generative_mem_mat'] = np.genfromtxt(
                path + '/generative_mem_mat.csv', delimiter=','
                ).astype(int)

        if os.path.isfile(path + '/clean_matrix.csv'):
            true_model_dict['clean_matrix'] = np.genfromtxt(
                path + '/clean_matrix.csv', delimiter=','
                ).astype(int)

        if os.path.isfile(path + '/parameters.csv'):
            parameters = np.genfromtxt(
                path + '/parameters.csv', delimiter=','
                )
            if parameters.size == 1:
                true_model_dict['number_of_errors'] = parameters[0]
            else:
                true_model_dict['generative_fbp'] = parameters[0]
                true_model_dict['generative_bgp'] = parameters[1]

        if n_patients is None:
            if len(true_model_dict) == 0:
                return(Dataset(matrix))
            else:
                return(Dataset(matrix, true_model_dict))
        else:
            matrix = matrix[:n_patients, :]
            if len(true_model_dict) == 0:
                return(Dataset(matrix))
            else:
                patients_stages = true_model_dict['patients_stages'][:n_patients]
                true_model_dict['patients_stages'] = patients_stages
                clean_matrix = true_model_dict['clean_matrix'][:n_patients, :]
                true_model_dict['clean_matrix'] = clean_matrix
                return(Dataset(matrix, true_model_dict))

    @classmethod
    def from_parameters(cls,
                        n_genes=25,
                        model_length=5,
                        n_genes_in_pws=None,
                        generative_fbp=0.05,
                        generative_bgp=0.05,
                        number_of_errors=None,
                        n_patients=1000,
                        number_of_passenger_genes=None,
                        PE=False,
                        debugging_mode=False,
                        seedno=None):

        if seedno is not None:
            seedno = np.random.seed(seedno)

        if n_genes_in_pws is None:
            if number_of_passenger_genes is None:
                splitters = np.random.choice(
                    n_genes-1,
                    model_length-1,
                    replace=False
                    )
                splitters += 1
                splitters.sort()
            else:
                splitters = np.random.choice(
                    n_genes-number_of_passenger_genes-1,
                    model_length-2,
                    replace=False
                    )
                splitters += 1
                splitters.sort()
                splitters = np.append(splitters, n_genes-number_of_passenger_genes)
        else:
            splitters = np.cumsum(n_genes_in_pws)

        if PE:
            stages_dist = np.ones(model_length)/(model_length - 1)
            stages_dist[-1] = 0
        else:
            stages_dist = np.ones(model_length)/model_length
        # Generating synthetic data:
        pathways = []
        gene_indices = np.arange(n_genes)
        if not debugging_mode:
            np.random.shuffle(gene_indices)

        pathways.append(gene_indices[0:splitters[0]])
        for i in np.arange(1, model_length-1):
            pathways.append(gene_indices[splitters[i-1]:splitters[i]])
        pathways.append(gene_indices[splitters[-1]:])

        generative_mem_mat = np.zeros((n_genes, model_length), dtype=np.int)
        for i in np.arange(model_length):
            for gene in pathways[i]:
                generative_mem_mat[gene, i] = 1

        patients_stages = np.random.choice(
            model_length,
            n_patients,
            p=stages_dist
            )
        dataset = np.zeros((n_patients, n_genes), dtype=np.int)
        for i in np.arange(n_patients):
            for j in np.arange(patients_stages[i]+1):
                gene = np.random.choice(pathways[j])
                dataset[i, gene] = 1
        clean_dataset = copy.deepcopy(dataset)
        true_model_dict = {
            'generative_mem_mat': generative_mem_mat,
            'patients_stages': patients_stages,
            'clean_matrix': clean_dataset
            }
        # Injecting the errors:
        if number_of_errors is None:
            for i in np.arange(n_patients):
                for gene in np.arange(n_genes):
                    if dataset[i, gene] == 1:
                        dataset[i, gene] = np.random.binomial(
                            1,
                            1 - generative_fbp
                            )
                    else:
                        dataset[i, gene] = np.random.binomial(
                            1,
                            generative_bgp
                            )
            true_model_dict['generative_fbp'] = generative_fbp
            true_model_dict['generative_bgp'] = generative_bgp
        else:
            errorous_bits = np.random.choice(
                n_patients*n_genes,
                number_of_errors,
                replace=False
                )
            for bit in errorous_bits:
                patient, gene = divmod(bit, n_genes)
                dataset[patient, gene] = 1 - dataset[patient, gene]
            true_model_dict['number_of_errors'] = number_of_errors
        return(cls(dataset, true_model_dict))

    @staticmethod
    def test():
        generative_fbp = 0.1
        generative_bgp = 0.01
        number_of_errors = 50000
        n_patients = 10000
        # First test:
        PE = False
        dataset = Dataset.from_parameters(
            generative_fbp=generative_fbp,
            generative_bgp=generative_bgp,
            n_patients=n_patients,
            PE=PE)
        diffs = dataset.clean_matrix-dataset.matrix
        n_flipbacks = np.sum(diffs > 0)
        n_clean_ones = np.sum(dataset.clean_matrix)
        fbp = n_flipbacks/n_clean_ones
        n_backgrounds = np.sum(diffs < 0)
        n_clean_zeros = n_patients*25 - n_clean_ones
        bgp = n_backgrounds/n_clean_zeros
        print(
            '\n test #1: passengers DO NOT EXIST, fixed ERROR RATE'
            )
        print(
            'generative fbp was {}, and the resulting fbp is {}'.format(
                generative_fbp,
                fbp
                )
            )
        print(
            'generative bgp was {}, and the resulting bgp is {}'.format(
                generative_bgp,
                bgp
                )
            )
        # Second test:
        PE = True
        dataset = Dataset.from_parameters(
            generative_fbp=generative_fbp,
            generative_bgp=generative_bgp,
            n_patients=n_patients,
            PE=PE)
        diffs = dataset.clean_matrix-dataset.matrix
        n_flipbacks = np.sum(diffs > 0)
        n_clean_ones = np.sum(dataset.clean_matrix)
        fbp = n_flipbacks/n_clean_ones
        n_backgrounds = np.sum(diffs < 0)
        n_clean_zeros = n_patients*25 - n_clean_ones
        bgp = n_backgrounds/n_clean_zeros
        print(
            '\n test #2: passengers EXIST, fixed ERROR RATE'
            )
        print(
            'generative fbp was {}, and the resulting fbp is {}'.format(
                generative_fbp,
                fbp
                )
            )
        print(
            'generative bgp was {}, and the resulting bgp is {}'.format(
                generative_bgp,
                bgp
                )
            )
        # Third test:
        PE = False
        dataset = Dataset.from_parameters(
            number_of_errors=number_of_errors,
            n_patients=n_patients,
            PE=PE)
        diffs = dataset.clean_matrix-dataset.matrix
        resulting_n_errors = np.sum(np.abs(diffs))
        print(
            '\n test #3: passengers DO NOT EXIST, fixed NUMBER OF ERRORS'
            )
        print(
            'number of errors was {}, and the resulting number is {}'.format(
                number_of_errors,
                resulting_n_errors
                )
            )
        # Third test:
        PE = True
        dataset = Dataset.from_parameters(
            number_of_errors=number_of_errors,
            n_patients=n_patients,
            PE=PE)
        diffs = dataset.clean_matrix-dataset.matrix
        resulting_n_errors = np.sum(np.abs(diffs))
        print(
            '\n Results of test #4: passengers EXIST, fixed NUMBER OF ERRORS'
                )
        print(
            'number of errors was {}, and the resulting number is {}'.format(
                number_of_errors,
                resulting_n_errors
                )
            )
        return()

    @staticmethod
    def calc_poco(mem_mat1, mem_mat2):
        n_genes = mem_mat1.shape[0]
        inp1 = np.zeros((n_genes, 1))
        inp2 = np.zeros_like(inp1)
        for i in range(0, n_genes):
            inp1[i] = np.where(mem_mat1[i, :] == 1)
            inp2[i] = np.where(mem_mat2[i, :] == 1)
        output = np.zeros((n_genes, n_genes))
        for first_gene in np.arange(n_genes-1):
            for second_gene in np.arange(first_gene+1, n_genes):
                if inp1[first_gene] > inp1[second_gene]:
                    if inp2[first_gene] > inp2[second_gene]:
                        output[first_gene, second_gene] = 1
                if inp1[first_gene] < inp1[second_gene]:
                    if inp2[first_gene] < inp2[second_gene]:
                        output[first_gene, second_gene] = 1
                if inp1[first_gene] == inp1[second_gene]:
                    if inp2[first_gene] == inp2[second_gene]:
                        output[first_gene, second_gene] = -1
        poco = np.sum(np.abs(output))/(n_genes*(n_genes-1)/2)
        return(poco)

    @staticmethod
    def calc_nocs(stages2, stages1):
        return (np.sum(stages1 == stages2))

if __name__ == "__main__":
    Dataset.test()
