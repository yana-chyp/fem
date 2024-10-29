import base_functions as base
import sympy as sp


def element_matrix(m):
    base_functions = base.get_base_functions(m)
    base_primes = [sp.diff(base_function, 'ksi') for base_function in base_functions]
    # print(base_primes)
    matrix = [[0 for j in range(m+1)] for i in range(m+1)]

    ksi = sp.symbols('ksi')
    for i in range(m+1):
        for j in range(m+1):
            value = sp.integrate(base_primes[i]*base_primes[j], (ksi, base.ksi_left, base.ksi_right))
            # print(i, j, value)
            matrix[i][j] = value



    return matrix