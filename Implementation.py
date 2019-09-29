from Functions import *

if __name__ == "__main__":
    initial_matrix = [[2, 5, 1], [3, 4, -2], [1, 1, 3]]
    free_terms = [[-2], [-3], [3]]

    print("\nInitial matrix:")
    print_matrix(initial_matrix, free_terms)

    # System triangularization
    transformed_matrix, transformed_free_terms = triangularize(initial_matrix, free_terms)
    print("\nTransformed matrix:")
    print_matrix(transformed_matrix, transformed_free_terms)

    # Searching for x vector
    x_vector = search_for_variables(transformed_matrix, transformed_free_terms)
    print("\nX vector: ", x_vector)

    # Searching for initial matrix determinant
    determinant = find_determinant(initial_matrix)
    print("\nDeterminant: ", determinant)

    # Searching for reverse matrix
    unit_matrix = create_unit_matrix(len(initial_matrix))
    transformed_matrix, transformed_unit_matrix = triangularize(initial_matrix, unit_matrix)
    reverse_matrix = search_for_variables(transformed_matrix, transformed_unit_matrix)
    print("\nReverse matrix:")
    print_matrix(reverse_matrix, None)

    # Verification
    print("\nVerification,  A * A^(-1) = \n", np.asarray(initial_matrix).dot(np.asarray(reverse_matrix)))

    # Incoherence
    print("\nIncoherence: \n", incoherence(initial_matrix, x_vector, free_terms))
