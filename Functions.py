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
    if not search_determinant:
        F = copy.deepcopy(free_terms)
    determinant = 1

    for i in range(0, len(A)):
        pivot_line_index = line_pivot_max(A, i)
        transposition(A, i, pivot_line_index)
        if not search_determinant:
            transposition(F, i, pivot_line_index)
        pivot = A[i][i]
        determinant *= pivot
        for j in range(i+1, len(A)):
            first = A[j][i]
            for k in range(0, len(A)):
                A[j][k] -= first * A[i][k] / pivot
            if not search_determinant:
                for k in range(0, len(F[0])):
                    F[j][k] -= first * F[i][k] / pivot

        for l in range(len(A)-1, -1, -1):
            A[i][l] /= pivot
        if not search_determinant:
            for l in range(len(F[0])-1, -1, -1):
                F[i][l] /= pivot

    if search_determinant:
        return determinant
    else:
        return A, F


def search_for_variables(A, F):
    n = len(A)
    X = [[0 for x in range(len(F[0]))] for i in range(n)]

    for i in range(0, len(F[0])):
        X[n-1][i] = F[n-1][i]
        for j in range(n-1, -1, -1):
            X[j][i] = F[j][i]
            for k in range(j+1, n):
                X[j][i] -= A[j][k]*X[k][i]

    return X


def find_determinant(matrix):
    return triangularize(matrix, search_determinant=True)


def print_matrix(matrix, free_terms_matrix = None):
    for row_number in range(len(matrix)):
        print("[", end="")
        for element in matrix[row_number]:
            print("% 9.5f " % element, end="")
        if free_terms_matrix:
            print("|", end="")
            for element in free_terms_matrix[row_number]:
                print(" % 9.5f" % element, end="")
        print("]")


def incoherence(matrix, solution_vector, free_terms):
    A = copy.deepcopy(matrix)
    X = copy.deepcopy(solution_vector)
    F = copy.deepcopy(free_terms)
    F_prime = np.asarray(A).dot(np.asarray(X))
    incoherence = [0] * len(F)
    for i in range(len(F)):
        incoherence[i] = (F[i] - F_prime[i]).tolist()
    return incoherence
