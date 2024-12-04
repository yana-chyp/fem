import cmath

import numpy as np
from matplotlib import pyplot as plt

import base_functions_2d as bs2d
import finite_element_2d as fe2d
import system_2d as s2d
import graph2d as g2d
import mesh_2d as m2d
def f(x, y):
    return 0

def exact_solution(x, y):
    return 31*x

def ug_1(x, y):
    return 0
def ug_2(x, y):
    return (-(y-0.5)**2 + 0.25)*8 + 2
def ug_3(x, y):
    return 2*x
def ug_4(x, y):
    return 2*x

b1 = 0
b2 = 0
d1 = 1
d2 = 1
p = 5
m = 5
degree = 3
element_type='D2QU4N'
ug = [ug_3, ug_2, ug_4, ug_1]

s2d.solve(b1, d1, b2, d2, p, m, degree, f, ug, element_type)

# base = bs2d.get_base_functions(degree)
# for f in base:
#     print(f)

# em = fe2d.element_matrix(degree)
# for row in em:
#     print('[' + ', '.join([f"{el:.3f}" for el in row]) + ']' )
#



# if degree == 1:
#     nodes, elements = m2d.uniform_mesh_level1(d1, d2, p, m, element_type)
# elif degree == 2:
#     nodes, elements = m2d.uniform_mesh_level2(d1, d2, p, m, element_type)
# elif degree == 3:
#     nodes, elements = m2d.uniform_mesh_level3(d1, d2, p, m, element_type)
# else:
#     print("unsupported degree")

# for i in range(len(nodes)): print(nodes[i])

# bounds = s2d.get_boundary_points(p, m)
# print(bounds)



# matrix = s2d.set_up_matrix(d1, d2, p, m, element_type, degree)
# print(50*"-")
# for row in matrix:
#     print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )

# print(50*"-")
# f_vec = s2d.set_up_vector(f, base, d1, d2, p, m, element_type, degree)
# print('[' + ', '.join([f"{el:.4f}" for el in f_vec]) + ']')

# matrix, f_vec = s2d.apply_boundary_conditions(matrix, f_vec, p, m, nodes, ug, degree)

# print(50*"-")
# print(50*"-")

# for row in matrix:
#     print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )
# print(50*"-")

# print('[' + ', '.join([f"{el:.4f}" for el in f_vec]) + ']')

# u = np.linalg.solve(matrix, f_vec)
# print('[' + ', '.join([f"{el:.4f}" for el in u]) + ']')

# m2d.plot_2d_solution(u, nodes, elements)
