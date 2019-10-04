from Functions import *

if __name__ == "__main__":
    k = 3
    alpha = 0.5 * k
    print("K =", k, ", alpha = ", alpha)

    initial_matrix = [[3.81, 0.25, 1.28, 0.75+alpha], [2.25, 1.32, 4.58 + alpha, 0.49],
                      [5.31, 6.28 + alpha, 0.98, 1.04], [9.39 + alpha, 2.45, 3.35, 2.28]]
    free_terms = [[4.21], [6.47 + alpha], [2.38], [10.48 + alpha]]

    print("\nInitial matrix:")
    print_matrix(initial_matrix, free_terms)

    # System triangularization
    transformed_matrix, transformed_free_terms, transformed_unit, determinant = triangularize(initial_matrix, free_terms)
    print("\nTransformed matrix:")
    print_matrix(transformed_matrix, transformed_free_terms)

    # Searching for x vector
    x_vector = search_for_variables(transformed_matrix, transformed_free_terms)
    print("\nX vector: ")
    print_matrix(x_vector)

    # Incoherence
    print("\nIncoherence: ")
    incoherence_matrix = incoherence(initial_matrix, x_vector, free_terms)
    for index in range(len(incoherence_matrix)):
        print("dx[%d] = " % index, incoherence_matrix[index][0])

    # Searching for initial matrix determinant
    print("\nDeterminant: ", determinant)

    # Searching for reverse matrix
    reverse_matrix = find_reverse_matrix(transformed_matrix, transformed_unit)
    print("\nReverse matrix:")
    print_matrix(reverse_matrix, None)

    # Verification
    print("\nVerification,  A * A^(-1) = \n", np.asarray(initial_matrix).dot(np.asarray(reverse_matrix)))
