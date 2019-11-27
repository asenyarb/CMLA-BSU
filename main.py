import time
from Implementation import *


if __name__ == "__main__":
    k = 3
    alpha = 0.5 * k
    initial = [[3.81, 0.25, 1.28, 0.75 + alpha], [2.25, 1.32, 4.58 + alpha, 0.49],
               [5.31, 6.28 + alpha, 0.98, 1.04], [9.39 + alpha, 2.45, 3.35, 2.28]]
    free_terms = [4.21, 6.47 + alpha, 2.38, 10.48 + alpha]
    print("Initial matrix:")
    Matrix.print(initial, free_terms)
    t_initial = [[7.6515, 2.2375, 2.262, 0.3675], [1.50, 7.53, -0.30, -1.21],
                 [2.25, 1.32, 6.08, 0.49], [3.81, 0.25, 1.28, 2.25]]
    t_free_terms = [8.4015, -1.83, 7.97, 4.21]

    iterative_methods = CMLAIterative(t_initial, t_free_terms, 10**5)
    print()
    t1 = time.perf_counter()
    x1, n1, B, b = iterative_methods.simple_iterative_method()
    t2 = time.perf_counter()
    print("Simple iterative method:\nx vector = ", x1, "\niterations number: ", n1)
    print("incoherence: ", iterative_methods.incoherence(x1), "B:", sep="\n")
    Matrix.print(B)
    print("b:", b)
    print("execution_time: ", t2 - t1, end="\n\n")
    t3 = time.perf_counter()
    x2, n2 = iterative_methods.jacobi_method()
    t4 = time.perf_counter()
    print("Jacobi method:\nTransformed matrix: ")
    Matrix.print(t_initial, t_free_terms)
    print("B:")
    Matrix.print(B)
    print("b:", b)
    print("||B||: ", Matrix.norm(B))
    print("||b||: ", Matrix.norm_for_vector(b))
    print("x vector = ", x2, "\niterations number: ", n2)
    print("incoherence: ", iterative_methods.incoherence(x2))
    print("execution_time: ", t4 - t3, end="\n\n")
    t5 = time.perf_counter()
    x4, n4 = iterative_methods.seidel_gauss_method()
    t6 = time.perf_counter()
    print("Seidel - Gauss method:\nx vector = ", x4, "\niterations number: ", n4)
    print("incoherence: ", iterative_methods.incoherence(x4))
    print("execution_time: ", t6-t5, end="\n\n")
