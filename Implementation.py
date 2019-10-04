from Functions import *

if __name__ == "__main__":
    k = 3
    alpha = 0.5 * k
    print("K = ", k, ", alpha = ", alpha)

    A = [[3.81, 0.25, 1.28, 0.75+alpha], [2.25, 1.32, 4.58 + alpha, 0.49],
                      [5.31, 6.28 + alpha, 0.98, 1.04], [9.39 + alpha, 2.45, 3.35, 2.28]]
    free_terms = [[4.21], [6.47 + alpha], [2.38], [10.48 + alpha]]

    print("\nLab_1 matrix:")
    print_matrix(A, free_terms)

    print("\nTransponated matrix:")
    transponated_matrix = transponate(A)
    print_matrix(transponated_matrix)

    print("\nInitial matrix:")
    initial_matrix = multiply_matrix(A, transponated_matrix)
    print_matrix(initial_matrix, free_terms)

    print("\nTransponated initial matrix * f")
    print_matrix(multiply_matrix(transponate(initial_matrix), free_terms))

    S, D, determinant = forward(initial_matrix, free_terms)
    print("\nB = D * S:")
    print_matrix(multiply_matrix(S, D))

  #  Y, X = backward(S, D, free_terms)
    print("\nY:")
   # print_matrix(Y)
    print("\nX:")
    #print_matrix(X)

    #print("\nDeterminant: ", determinant)
