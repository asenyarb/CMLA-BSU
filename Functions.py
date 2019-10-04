import copy
import numpy as np
from mpmath import *


def transponate(matrix):
    transponated_matrix = [[0 for i in range(len(matrix))] for i in range(len(matrix))]
    for i in range(len(transponated_matrix)):
        for j in range(len(transponated_matrix[0])):
            transponated_matrix[i][j] = matrix[j][i]
    return transponated_matrix


def multiply_matrix(A, B):
    mp.dps = 1000
    result = [[mpf("0") for i in range(len(B[0]))] for j in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] = fmul(A[i][k], B[k][j])
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
    D = [0] * len(initial_matrix)
    determinant = 1

    for i in range(len(S)):
        undersquare_value = initial_matrix[i][i]
        for k in range(i-1):
            a = fmul(S_prime[i][k], S[k][i])
            b = fmul(a, D[k])
            undersquare_value = fsub(undersquare_value, b)
        D[i] = sign(undersquare_value)
        S_prime[i][i] = S[i][i] = sqrt(fmul(undersquare_value, D[i]))
        determinant *= fmul(power(S[i][i], 2), D[i])
        for j in range(i+1, len(S[0])):
            nominator = initial_matrix[i][j]
            for k in range(i-1):
                nominator = fsub(nominator, fmul(S_prime[k][i], fmul(S[k][j], D[k])))
            S[i][j] = S_prime[j][i] = fdiv(nominator, fmul(D[i], S[i][i]))
    return S, D, determinant

def backward(S, D, free_terms):
    mp.dps = 1000
    X = [[0] for i in range(len(S))]
    Y = [[0] for i in range(len(S))]

    # Find Y matrix
    for i in range(len(S)):
        nominator = free_terms[i][0]
        for k in range(i-1):
            nominator = fsub(nominator, fmul(S[k][i], Y[k][0]))
        Y[i][0] = fdiv(nominator, S[i][i])

    # Find X matrix
    for i in range(len(S)-1, -1, -1):
        nominator = Y[i][0]
        for k in range(j+1, len(S)):
            nominator = fsub(nominator, fmul(D[j], fmul(S[j][k], X[k])))
        X[i][0] = fdiv(nominator, fmul(D[j], S[j][j]))
    return Y, X