import sympy as sp
import numpy as np
import nastia as nastia
import matplotlib.pyplot as plt


def get_base_function(i, mesh, side):
    # def N_i(x) :
    x = sp.symbols('x')
    expression = 0
    if i == 0:
        # if sp.And(mesh[0] <= x, x <= mesh[1]):
        if side==2:
            expression = (mesh[1] - x) / (mesh[1] - mesh[0])
        else:
           expression = 0
    elif i == len(mesh) - 1:
         # if sp.And(mesh[-2] <= x, x <= mesh[-1]):
         if side==1:
            expression = (x - mesh[-2]) / (mesh[-1] - mesh[-2])
         else:
            expression = 0
    else:
        # if sp.And(mesh[i-1] <= x, x <= mesh[i]):
        if side==1:             #let this mean left side, elif - right
            expression =  (x - mesh[i-1]) / (mesh[i] - mesh[i-1])
        # elif sp.And(mesh[i] <= x, x <= mesh[i+1]):
        elif side==2:
            expression = (mesh[i+1] - x) / (mesh[i+1] - mesh[i])
        else:
            expression = 0
    return expression
    # func = sp.lambdify(x, expression)
    # return func

def set_up_element_matrix(i, mesh) :
    N_1 = get_base_function(i, mesh, 2)
    N_2 = get_base_function(i+1, mesh, 1)

    N_1_prime = sp.diff(N_1, 'x')
    N_2_prime = sp.diff(N_2, 'x')

    # x = sp.symbols('x');
    x_1 = mesh[i]; x_2 = mesh[i+1]
    x = sp.symbols('x')
    matrix_11 = sp.integrate(N_1_prime*N_1_prime, (x, x_1, x_2));
    matrix_12 = sp.integrate(N_1_prime*N_2_prime, (x, x_1, x_2));
    matrix_21 = sp.integrate(N_2_prime*N_1_prime, (x, x_1, x_2));
    matrix_22 = sp.integrate(N_2_prime*N_2_prime, (x, x_1, x_2));
    matrix = np.array([np.array([matrix_11, matrix_12]),
                    np.array([matrix_21, matrix_22])])
    return matrix

def set_up_general_matrix(mesh):
    matrix = []
    row = [0 for i in range(len(mesh))]
    for i in range(len(mesh)-1):
        # print('i = ', i)
        K_i = set_up_element_matrix(i, mesh)
        # print("K_i = ", [[f"{el:.2f}" for el in K_i_row] for K_i_row in K_i])
        # np.set_printoptions(precision=2, suppress=True)
        # print(K_i)
        for j in range(len(K_i)):
            row[i+j] += K_i[0][j]
        # print("row = ", [f"{el:.2f}" for el in row])
        matrix.append(row[:])
        # print("matrix = ", [[f"{el:.2f}" for el in mat_row] for mat_row in matrix])

        # row = [0 for i in mesh.len() + 1]
        if i>0:
            row[i-1] = 0
        for j in range(len(K_i)):
            row[i+j] = K_i[1][j]
        # print("row = ", row)
        # print("matrix = ", matrix)
        # print("row = ", [f"{el:.2f}" for el in row])
        # print("matrix = ", [[f"{el:.2f}" for el in mat_row] for mat_row in matrix])

    matrix.append(row[:])
    K = np.array(matrix, dtype=np.float64)
    return K;

def use_boundary_conditions(K, f_arr, c, d):
    n = len(K)
    # red_K = [[0 for _ in range(n-2)] for _ in range(n-2)]
    # for i in range(n-2): red_K.append([])
    # for j in range(n-2):
        # red_K[0][j] = (K[0][j+1] + K[1][j+1])
    # for i in range(1, n-3):
    #     for j in range(n-2):
    #         red_K[i][j] = (K[i+1][j+1])
    # for j in range(n-2):
    #     red_K[n-3][j] = (K[n-2][j]+K[n-1][j])

    # for i in  range(n):
    #     f_arr[i] -= K[i][0]*c - K[i][n-1]*d

    # f_arr[1] += f_arr[0]; f_arr[0] = 0
    # f_arr[n-2] += f_arr[n-1]; f_arr[n-1] = 0
        # if i>0 & i<n-1 & i%3==0:
        #     f_arr[i]*=2
    # return red_K;

    red_K = [[0 for _ in range(n)] for _ in range(n)]
    red_K[0][0] = 1; f_arr[0] = c;
    red_K[n-1][n-1] = 1; f_arr[n-1] = d;
    for i in range(1, n-1):
        for j in range(1, n-1):
            red_K[i][j] = (K[i][j])

    return np.array(red_K, dtype=np.float64)

def get_f_arr(mesh):
    f_arr = []
    for i in range(len(mesh)):
        # integral_phi = nastia.integrate_phi_i(i, mesh)
        integral_f_phi = nastia.integrate_f_phi_i(i, mesh)
        f_arr.append(integral_f_phi)
    return np.array(f_arr, dtype=np.float64)
    # return f_arr

def solve_system(red_K, f_arr):
    u = np.linalg.solve(red_K, f_arr)
    return u;

def plot_solution(u, mesh, a, b):
    x_values = np.linspace(a, b, 100)
    approx_values = [nastia.approximate_solution(x, mesh, u) for x in x_values]
    exact_values = [nastia.exact_solution(x) for x in x_values]
    plt.plot(x_values, approx_values, label='Наближений розв\'язок')
    plt.plot(x_values, exact_values, label='Точний розв\'язок', linestyle='--')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('Наближений розв\'язок на інтервалі [0, 1]')
    plt.legend()
    plt.grid(True)
    plt.show()

# def solve_by_fem()