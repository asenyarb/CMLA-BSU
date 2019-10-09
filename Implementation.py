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
    transposed_matrix = transponate(A)
    print_matrix(transposed_matrix)

    # Multiply initial matrix by transposed initial matrix to get symmetrical matrix
    print("\nInitial matrix:")
    initial_matrix = multiply_matrix(transposed_matrix, A)
    print_matrix(initial_matrix)

    # We also should multiply free terms matrix if we want our new equation to be equivalent to the initial one
    print("\nTransponated initial matrix * f")
    free_terms = multiply_matrix(transposed_matrix, free_terms)
    print_matrix(free_terms)

    S, D, determinant, Y = forward(initial_matrix, free_terms)
    print("\nB = D * S:")
    print_matrix(multiply_matrix(D, S))
    print("\nY:")
    print_matrix(Y)

    X = backward(S, D, Y)
    print("\nX:")
    print_matrix(X)

    print("\n|det(A)| = %.10f" % determinant)
