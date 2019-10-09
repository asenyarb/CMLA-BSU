import copy
import numpy as np
from mpmath import *


def transponate(matrix):
    transposed_matrix = [[0 for i in range(len(matrix))] for i in range(len(matrix))]
    for i in range(len(transposed_matrix)):
        for j in range(len(transposed_matrix[0])):
            transposed_matrix[i][j] = matrix[j][i]
    return transposed_matrix


def multiply_matrix(A, B):
    mp.dps = 1000
    result = [[mpf("0") for i in range(len(B[0]))] for j in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            result[i][j] = mpf("0")
            for k in range(len(B)):
                result[i][j] += fmul(A[i][k], B[k][j])
    return result


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


def forward(initial_matrix, free_terms):
    mp.dps = 1000
    S = [[mpf("0") for i in range(len(initial_matrix[0]))] for i
         in range(len(initial_matrix))]
    S_prime = [[mpf("0") for i in range(len(initial_matrix[0]))] for i
               in range(len(initial_matrix))]
    D = [[mpf("0") for i in range(len(initial_matrix[0]))] for j in range(len(initial_matrix))]
    determinant = 1

    for i in range(len(S)):
        undersquare_value = initial_matrix[i][i]
        for k in range(i):
            undersquare_value = fsub(undersquare_value, fmul(power(S[k][i], 2), D[k][k]))
        D[i][i] = sign(undersquare_value)
        S_prime[i][i] = S[i][i] = sqrt(fmul(undersquare_value, D[i][i]))
        determinant *= fmul(power(S[i][i], 2), D[i][i])
        divisor = fmul(D[i][i], S[i][i])
        for j in range(i+1, len(S[0])):
            nominator = initial_matrix[i][j]
            for k in range(i):
                nominator = fsub(nominator, fmul(S[k][i], fmul(S[k][j], D[k][k])))
            S[i][j] = S_prime[j][i] = fdiv(nominator, divisor)

    Y = [[mpf("0")] for i in range(len(S))]
    # Find Y matrix
    for i in range(len(S)):
        nominator = free_terms[i][0]
        for k in range(i):
            nominator = fsub(nominator, fmul(S[k][i], Y[k][0]))
        Y[i][0] = fdiv(nominator, S[i][i])

    return S, D, sqrt(determinant), Y

def backward(S, D, Y):
    mp.dps = 1000
    X = [[mpf("0")] for i in range(len(S))]
    # Find X matrix
    for j in range(len(S)-1, -1, -1):
        nominator = Y[j][0]
        for k in range(j+1, len(S)):
            nominator = fsub(nominator, fmul(D[j][j], fmul(S[j][k], X[k][0])))
        X[j][0] = fdiv(nominator, fmul(D[j][j], S[j][j]))
    return X