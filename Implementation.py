from copy import deepcopy
import sys


class Matrix:
    @staticmethod
    def norm(matrix):
        if matrix:
            norm = sys.float_info.min
            for i in range(len(matrix)):
                somme = 0
                for j in range(len(matrix[0])):
                    somme += abs(matrix[i][j])
                if somme > norm:
                    norm = somme
            return norm
        else:
            return None

    @staticmethod
    def norm_for_vector(vector):
        norm = 0
        for el in vector:
            norm += el
        return norm

    @staticmethod
    def multiply(matrix1, matrix2):
        result = [[0 for _ in range(len(matrix2[0]))] for _ in range(len(matrix1))]
        for i in range(len(matrix1)):
            for j in range(len(matrix2[0])):
                result[i][j] = 0
                for k in range(len(matrix2)):
                    result[i][j] += matrix1[i][k] * matrix2[k][j]
        return result

    @staticmethod
    def multiply_matrix_column(matrix, column):
        result = [0] * len(column)
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                result[i] += matrix[i][j] * column[j]
        return result

    @staticmethod
    def sub(m1, m2):
        return [[m1[i][j] - m2[i][j]
                 for j in range(len(m1[0]))]
                for i in range(len(m1))]

    @staticmethod
    def add(m1, m2):
        return [[m1[i][j] + m2[i][j]
                 for j in range(len(m1[0]))]
                for i in range(len(m1))]

    @staticmethod
    def add_v(v1, v2):
        return [v1[i] + v2[i] for i in range(len(v1))]

    @staticmethod
    def div(m, a):
        matr = [[m[i][j] / a
                 for j in range(len(m[0]))]
                for i in range(len(m))]
        return matr

    @staticmethod
    def transposed(matrix):
        return [[matrix[j][i] for j in range(len(matrix[0]))]
                for i in range(len(matrix))]

    @staticmethod
    def print(matrix, free_terms=None):
        free_terms_str = ["" for _ in range(len(matrix))]\
            if not free_terms else ["| " + str(f_t) + " " for f_t in free_terms]
        for i in range(len(matrix)):
            print("[ " + " ".join(list(map(str, matrix[i]))) + free_terms_str[i] + "]")


def max_difference(x_k1, x_k):
    diff_list = [abs(x_k1[i] - x_k[i]) for i in range(len(x_k1))]
    return max(diff_list)


class CMLAIterative:
    initial_matrix = []
    lines_number = 0
    columns_number = 0
    precision = 0
    free_terms = []
    unit_matrix = []

    def __init__(self, matrix, free_terms, precision):
        self.initial_matrix = deepcopy(matrix)
        self.free_terms = deepcopy(free_terms)
        self.lines_number = len(matrix)
        self.columns_number = 0 if not self.lines_number else len(matrix[0])
        self.unit_matrix = [[0 if i != j else 1 for j in range(len(self.initial_matrix[0]))]
                            for i in range(len(self.initial_matrix))]
        if precision:
            self.precision = 1 / precision
        else:
            raise EnvironmentError("You have to specify the precision value > 0!")

    def set_precision(self, precision):
        self.precision = 1 / precision

    @staticmethod
    def util_sum(matrix, i, x):
        s = 0
        for j in range(len(matrix[0])):
            if j != i:
                s += matrix[i][j] * x[j]
        return s

    def simple_iterative_method(self):
        # A transposed * A
        tr = Matrix.transposed(self.initial_matrix)
        matrix_a_t_a = Matrix.multiply(tr, self.initial_matrix)
        a_t_a_norm = Matrix.norm(matrix_a_t_a)
        matrix_b = Matrix.sub(self.unit_matrix, Matrix.div(matrix_a_t_a, a_t_a_norm))
        vector_b = Matrix.multiply_matrix_column(tr, [self.free_terms[i] / a_t_a_norm for i in range(len(self.free_terms))])
        iterations_number = 0
        x_prev = [0 for _ in range(len(self.free_terms))]
        x_vector = deepcopy(vector_b)
        while max_difference(x_vector, x_prev) > self.precision:
            x_prev = deepcopy(x_vector)
            x_vector = Matrix.add_v(Matrix.multiply_matrix_column(matrix_b, x_vector), vector_b)
            iterations_number += 1
        return x_vector, iterations_number, matrix_b, vector_b

    def jacobi_method(self):
        iterations_number = 0
        x_prev = [0 for _ in range(len(self.free_terms))]
        x_vector = deepcopy(self.free_terms)
        while max_difference(x_vector, x_prev) > self.precision:
            x_prev = deepcopy(x_vector)
            x_vector = [(-self.util_sum(self.initial_matrix, i, x_vector) +
                         self.free_terms[i]) / self.initial_matrix[i][i]
                        for i in range(len(x_vector))]
            iterations_number += 1
        return x_vector, iterations_number

    @staticmethod
    def seidel_sum(b_matrix, i, x):
        s = 0
        for j in range(len(x)):
            if i != j:
                s += b_matrix[i][j] * x[j]
        return s

    def seidel_gauss_method(self):
        x_prev = [0 for _ in range(len(self.free_terms))]
        x_vector = deepcopy(self.free_terms)
        iterations_number = 0
        while max_difference(x_vector, x_prev) > self.precision:
            x_prev = deepcopy(x_vector)
            for i in range(len(x_vector)):
                x_vector[i] = (-self.seidel_sum(self.initial_matrix, i, x_vector) +
                               self.free_terms[i]) / self.initial_matrix[i][i]
            iterations_number += 1
        return x_vector, iterations_number

    def incoherence(self, x_vector):
        result = Matrix.multiply_matrix_column(self.initial_matrix, x_vector)
        return [self.free_terms[i] - result[i]
                for i in range(len(self.free_terms))]
