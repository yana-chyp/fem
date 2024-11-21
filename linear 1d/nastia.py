import cmath
import sympy as sp
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

#division into segments
def get_mesh(a, b, n):
    h = (b - a) / n
    arr = [a + i * h for i in range(n + 1)]
    return arr
#creation of basis functions (for subsequent calculations of F)
def phi_i(x, i, mesh):
    if i == 0:
        if mesh[0] <= x <= mesh[1]:
            return (mesh[1] - x) / (mesh[1] - mesh[0])
        else:
            return 0
    elif i == len(mesh) - 1:
        if mesh[-2] <= x <= mesh[-1]:
            return (x - mesh[-2]) / (mesh[-1] - mesh[-2])
        else:
            return 0
    else:
        if mesh[i-1] <= x <= mesh[i]:
            return (x - mesh[i-1]) / (mesh[i] - mesh[i-1])
        elif mesh[i] <= x <= mesh[i+1]:
            return (mesh[i+1] - x) / (mesh[i+1] - mesh[i])
        else:
            return 0
#no function needed, was to check phi
def integrate_phi_i(i, mesh):
    if i == 0:
        integral, error = quad(lambda x: phi_i(x, i, mesh), mesh[0], mesh[1])
    elif i == len(mesh) - 1:
        integral, error = quad(lambda x: phi_i(x, i, mesh), mesh[-2], mesh[-1])
    else:
        integral_left, error_left = quad(lambda x: phi_i(x, i, mesh), mesh[i-1], mesh[i])
        integral_right, error_right = quad(lambda x: phi_i(x, i, mesh), mesh[i], mesh[i+1])
        integral = integral_left + integral_right
    return integral
#Vector F
def integrate_f_phi_i(i, mesh):
    if i == 0:
        integral, error = quad(lambda x: f(x) * phi_i(x, i, mesh), mesh[0], mesh[1])
    elif i == len(mesh) - 1:
        integral, error = quad(lambda x: f(x) * phi_i(x, i, mesh), mesh[-2], mesh[-1])
    else:
        integral_left, error_left = quad(lambda x: f(x) * phi_i(x, i, mesh), mesh[i-1], mesh[i])
        integral_right, error_right = quad(lambda x: f(x) * phi_i(x, i, mesh), mesh[i], mesh[i+1])
        integral = integral_left + integral_right
    return integral
#approximate solution
def approximate_solution(x, mesh, u):
    result = 0
    for i in range(len(mesh)):
        result += u[i] * phi_i(x, i, mesh)
    return result
#function
def f(x):
    # return 1
    # return -sp.exp(x)
    # return -12*x**2
    return sp.sin(x)

#solution -u'' =1 original function
def exact_solution(x):
    # return 0.5 * (x - x**2)
    # return sp.exp(x) + x
    # return x**4 - 5*x + 15
    return sp.sin(x) + x

# a = 0
# b = 1
# n = 2
#
# mesh = get_mesh(a, b, n)
# print("Вузли: ", mesh)
# K = np.array([[-1, -11, -3],
#               [1, 1, 0],
#               [2, 5, 1]])
#
# f_arr = []
# for i in range(len(mesh)):
#     integral_phi = integrate_phi_i(i, mesh)
#     print(f"Інтеграл від phi_{i} : {integral_phi:.6f}")
#
#     integral_f_phi = integrate_f_phi_i(i, mesh)
#     f_arr.append(integral_f_phi)
#     print(f"Інтеграл f(x) * phi_{i}(x) : {integral_f_phi:.6f}")
#
# print("\nМасив f_arr з інтегралами:")
# print(f_arr)
# #no function needed, was to check the result
# F = np.array([-37, -1, 10])
# #the result of finding the unknown vector (K (stiffness matrix)*u =F )
# u = np.linalg.solve(K, f_arr)
# u2 = np.linalg.solve(K, F)
#
# print("Рішення u:", u)
# print("Рішення u2:", u2)
# #we look for many points to build a picture
# x_values = np.linspace(a, b, 100)
# approx_values = [approximate_solution(x, mesh, u2) for x in x_values]
# exact_values = [exact_solution(x) for x in x_values]
# # plt.plot(x_values, approx_values, label='Наближений розв\'язок')
# # plt.xlabel('x')
# # plt.ylabel('u(x)')
# # plt.title('Наближений розв\'язок на інтервалі [0, 1]')
# # plt.legend()
# # plt.grid(True)
# # plt.show()
#
#
# #draw a graph
# plt.plot(x_values, approx_values, label='Наближений розв\'язок')
# plt.plot(x_values, exact_values, label='Точний розв\'язок', linestyle='--')
# plt.xlabel('x')
# plt.ylabel('u(x)')
# plt.title('Наближений і точний розв\'язки на інтервалі [0, 1]')
# plt.legend()
# plt.grid(True)
# plt.show()
