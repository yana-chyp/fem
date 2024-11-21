from matplotlib import pyplot as plt

import finite_element as finel
import base_functions as base
import sympy as sp
import nastia as nastia

def set_up_matrix(m, mesh):
    n = len(mesh) - 1
    dim = n*m+1
    matrix = [[0 for j in range(dim)] for i in range(dim)]
    # matrix = []
    element = finel.element_matrix(m)

    for i in range(n):
        J = (mesh[i+1]-mesh[i])/(base.ksi_right-base.ksi_left)
        for k in range(m+1):
            for l in range(m+1):
                matrix[i*m+k][i*m+l] += element[k][l] / J
    return matrix

def apply_boundary_conditions(matrix, f_arr, fa, fb):
    dim = len(matrix)
    matrix[0][0] = 1
    #could be to m, not to dim
    for j in range(1, dim): matrix[0][j] = 0
    matrix[-1][-1] = 1
    for j in range(dim-1): matrix[-1][j] = 0
    f_arr[0] = fa
    f_arr[-1] = fb

def set_up_vector(f, base_functions, mesh):
    n = len(mesh) - 1
    m = len(base_functions) - 1
    dim = n*m+1
    f_arr = [0 for i in range(dim)]
    ksi = sp.symbols('ksi')
    x = sp.symbols('x')

    for i in range(n):
        # f = f(base.get_x(ksi, i, mesh))
        # lin = sp.sympify(base.get_x(ksi, i ,mesh))
        lin_x = (ksi+1)/2 * (mesh[i+1] - mesh[i]) + mesh[i]
        lin_ksi = 2*(x - mesh[i])/(mesh[i+1] - mesh[i]) - 1
        f_expr = sp.sympify(nastia.f(x))
        # f_expr = f_expr.subs(x, lin_x)
        J = (mesh[i+1]-mesh[i])/(base.ksi_right-base.ksi_left)
        for k in range(m+1):
            phi = base_functions[k].subs('ksi', lin_ksi)
            # f_arr[i*m+k] = sp.integrate(f_expr * base_functions[k], (ksi, base.ksi_left, base.ksi_right)) * J
            f_arr[i*m+k] += sp.integrate(f_expr * phi, (x, mesh[i], mesh[i+1]))

    return f_arr

def get_solution(q, base_functions, mesh):
    x = sp.symbols('x')
    ksi = sp.symbols('ksi')
    m = len(base_functions) - 1
    n = len(mesh) - 1
    solution = 0
    for i in range(n):
        lin = 2*(x - mesh[i])/(mesh[i+1] - mesh[i]) - 1
        J = (mesh[i+1]-mesh[i])/(base.ksi_right-base.ksi_left)
        for k in range(m+1):
            solution += q[i*m+k]*base_functions[k].subs(ksi, lin)
    return solution

def get_interval(x, mesh):
    # index = -1
    for i in range(len(mesh)-1):
        if mesh[i] <= x <= mesh[i+1]: return i
    return -1
    # return index

def get_u_value(x, mesh, q, base_functions):
    m = len(base_functions)-1
    i = get_interval(x, mesh)
    solution = 0
    for k in range(m+1):
        func = sp.lambdify(sp.symbols('ksi'), base_functions[k])
        ksi = 2*(x - mesh[i])/(mesh[i+1] - mesh[i]) - 1
        # print("i = ", i, " ksi = ", x, " func(ksi) = ", func(ksi))
        solution += q[i*m+k]*func(ksi)
    return solution

def plot_solution(u_values, x_values):
    exact_values = [nastia.exact_solution(x) for x in x_values]
    plt.plot(x_values, exact_values, label='Точний розв\'язок', linestyle='--')
    plt.plot(x_values, u_values, label='Наближений розв\'язок')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('Наближений розв\'язок на інтервалі [0, 1]')
    plt.legend()
    plt.grid(True)
    plt.show()
