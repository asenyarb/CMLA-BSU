import copy
import numpy as np


def create_unit_matrix(dimension):
    unit_matrix = [[0 for i in range(dimension)] for i in range(dimension)]
    for i in range(dimension):
        unit_matrix[i][i] = 1
    return unit_matrix


def line_pivot_max(matrix, c):
    max_element = matrix[c][c]
    max_index = c
    for i in range(c, len(matrix)):
        if matrix[i][c] > max_element:
            max_element = matrix[i][c]
            max_index = i
    return max_index


def transposition(matrix, line_1, line_2):
    temporary = matrix[line_1]
    matrix[line_1] = matrix[line_2]
    matrix[line_2] = temporary


def triangularize(matrix, free_terms=None, search_determinant=False):
    A = copy.deepcopy(matrix)
    F = copy.deepcopy(free_terms)
    unit_matrix = create_unit_matrix(len(A))
    determinant = 1

    for i in range(0, len(A)):
        pivot_line_index = line_pivot_max(A, i)
        if i != pivot_line_index:
            transposition(A, i, pivot_line_index)
            transposition(F, i, pivot_line_index)
            transposition(unit_matrix, i, pivot_line_index)
            determinant *= -1
        pivot = A[i][i]
        determinant *= pivot
        for j in range(i+1, len(A)):
            first = A[j][i]
            for k in range(0, len(A)):
                A[j][k] -= first * A[i][k] / pivot
            F[j] -= first * F[i] / pivot
            for k in range(0, len(unit_matrix[0])):
                unit_matrix[j][k] -= first * unit_matrix[i][k] / pivot

        for l in range(len(A)-1, -1, -1):
            A[i][l] /= pivot
        F[i] /= pivot
        for l in range(len(unit_matrix[0]) - 1, -1, -1):
            unit_matrix[i][l] /= pivot

    return A, F, unit_matrix, determinant


def search_for_variables(A, F):
    n = len(A)
    X = [0 for i in range(n)]

    X[n-1] = F[n-1]
    for j in range(n-1, -1, -1):
        X[j] = F[j]
        for k in range(j+1, n):
            X[j] -= A[j][k]*X[k]

    return X


def print_matrix(matrix, free_terms_matrix=None, free_terms_vector=None):
    for row_number in range(len(matrix)):
        print("[", end="")
        for element in matrix[row_number]:
            print("% 9.5f " % element, end="")
        if free_terms_matrix:
            print("|", end="")
            for element in free_terms_matrix[row_number]:
                print(" % 9.5f" % element, end="")
        elif free_terms_vector:
            print("|", end="")
            print(" % 9.5f" % free_terms_vector[row_number], end="")
        print("]")


def print_vector(vector):
    for el in vector:
        print("[ %8.5f ]" % el)


def incoherence(matrix, solution_vector, free_terms):
    A = copy.deepcopy(matrix)
    X = copy.deepcopy(solution_vector)
    F = copy.deepcopy(free_terms)
    F_prime = np.asarray(A).dot(np.asarray(X))
    incoherence = [0] * len(F)
    for i in range(len(F)):
        incoherence[i] = (F[i] - F_prime[i]).tolist()
    return incoherence


def find_reverse_matrix(triangularized_A, triangularized_unit):
    reverse_matrix = [[0 for x in range(len(triangularized_unit[0]))] for i
                      in range(len(triangularized_unit))]
    for i in range(len(triangularized_unit)):
        column_unit = [triangularized_unit[j][i] for j in range(len(triangularized_unit))]
        x_column = search_for_variables(triangularized_A, column_unit)
        for j in range(len(reverse_matrix)):
            reverse_matrix[j][i] = x_column[j]
    return reverse_matrix