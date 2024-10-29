import numpy as np

import nastia as nastia
import fem_lib as fem
import sympy as sp

from fem_lib import use_boundary_conditions

a = 0
b = 1
n = 10
c = nastia.exact_solution(a)
d = nastia.exact_solution(b)


mesh = nastia.get_mesh(a, b, n)
print("Вузли: ", mesh)


# K_1 = fem.set_up_element_matrix(0, mesh)
# sp.pprint(K_1)

K = fem.set_up_general_matrix(mesh)
print("K: ")
# np.set_printoptions(precision=2)
# sp.pprint(K)
for row in K:
    print('[' + ', '.join([f"{el:.2f}" for el in row]) + ']' )

# print(K)

f_arr = fem.get_f_arr(mesh)
print("f_arr = ", f_arr)
red_K = use_boundary_conditions(K, f_arr, c, d)
print("red_K = ")
for row in red_K:
    print('[' + ', '.join([f"{el:.2f}" for el in row]) + ']' )

u = fem.solve_system(red_K, f_arr)
# u = np.insert(u, 0, c)
# u = np.append(u, d)


print("mesh = ", mesh)
print("u(x_i) = ", u)

# print("difference: ", K*u-f_arr);

fem.plot_solution(u, mesh, a, b)