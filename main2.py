import base_functions as base
import finite_element as finel
import nastia as nastia
from fem_lib import get_base_function
import matplotlib.pyplot as plt
import sympy as sp
import system_of_elements as system
import fem_lib as fem


# x = sp.symbols('x')


a = -1
b = 1
c = nastia.exact_solution(a)
d = nastia.exact_solution(b)
n = 3
m = 1
mesh = nastia.get_mesh(a, b, n)
base_functions = base.get_base_functions(m)
matrix = system.set_up_matrix(m, mesh)

for row in matrix:
    print('[' + ', '.join([f"{el:.2f}" for el in row]) + ']' )
print("-"*30)

f_arr = system.set_up_vector(nastia.f, base_functions, mesh)
system.apply_boundary_conditions(matrix, f_arr, c, d)

for row in matrix:
    print('[' + ', '.join([f"{el:.2f}" for el in row]) + ']' )
print("-"*30)

print('[' + ', '.join([f"{el:.2f}" for el in f_arr]) + ']' )
print("-"*30)


q = fem.solve_system(matrix, f_arr)
print("q(x_i) = ", q)
print("-"*30)
print("-"*30)
print("-"*30)


# x = sp.symbols('x')
# for i in range(n):
#     for k in range(m):
#         print("i = ", i, "mesh[i] = ", mesh[i], "k = ", k)
#         print(sp.lambdify(x, base_functions[k].subs('ksi', 2*(x - mesh[i])/(mesh[i+1] - mesh[i]) - 1))(mesh[i]) )
#     print('\n')

u = system.get_solution(q, base_functions, mesh)
print("u(x) = ", u)
print("-"*30)

u = sp.lambdify(sp.symbols('x'), u)
u_values = [u(x) for x in mesh]

print("mesh = ", mesh)
print("exact solution = ", [nastia.exact_solution(x) for x in mesh])

# fem.plot_solution(u_values, mesh, a, b)



# ksi_element = finel.element_matrix(m)
# for i in range(m + 1):
#     print(ksi_element[i])

# div = 100
#
# h = (base.ksi_right - base.ksi_left) / div
# ksi_s = [base.ksi_left + k*h for k in range(div+1)]
# # print(ksi_s)
# for i in range(n) :
#
#     j = 0
#     for expression in base_functions:
#         # print(expression)
#         function = sp.lambdify('ksi', expression)
#         x = [base.get_x(ksi, i, mesh) for ksi in ksi_s]
#         y = [function(ksi) for ksi in ksi_s]
#         plt.plot(x, y, label="phi[%i][%i]"%(i, j))
#         j+=1
#
# plt.xlabel('x')
# plt.ylabel('y')
# plt.legend()
# plt.grid(True)
# plt.show()