import numpy as np
from matplotlib import pyplot as plt

import base_functions_2d as bs2d
import finite_element_2d as fe2d
import system_2d as s2d
import graph2d as g2d

def f(x, y):
    return 0

def exact_solution(x, y):
    return 10*x

def ug_1(x, y):
    return 0
def ug_2(x, y):
    return 10
def ug_3(x, y):
    return 10*x
def ug_4(x, y):
    return 10*x


base = bs2d.get_base_functions()
# for f in base_functions:
#     print(f)

# em = fe2d.element_matrix()
# for row in em:
#     # print(row)
#     print('[' + ', '.join([f"{el:.3f}" for el in row]) + ']' )
#

d1 = 1
d2 = 1
p = 14
m = 13
element_type='D2QU4N'

nodes, elements = g2d.uniform_mesh(d1, d2, p, m, element_type)

# bounds = s2d.get_boundary_points(p, m)
# print(bounds)



matrix = s2d.set_up_matrix(d1, d2, p, m)
print(50*"-")
for row in matrix:
    print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )

print(50*"-")
f_vec = s2d.set_up_vector(f, base, d1, d2, p, m, element_type)
print('[' + ', '.join([f"{el:.4f}" for el in f_vec]) + ']')

ug = [ug_3, ug_2, ug_4, ug_1]
matrix, f_vec = s2d.apply_boundary_conditions(matrix, f_vec, p, m, nodes, ug)

print(50*"-")
print(50*"-")

for row in matrix:
    print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )
print(50*"-")

print('[' + ', '.join([f"{el:.4f}" for el in f_vec]) + ']')


u = np.linalg.solve(matrix, f_vec)
print('[' + ', '.join([f"{el:.4f}" for el in u]) + ']')


g2d.plot_2d_solution(u, nodes, elements)

# NL = nodes
# EL = elements
# NoN = NL.shape[0]
# NoE = EL.shape[0]

# plt.figure(figsize=(8, 8))
#
#     # Анотація вузлів
# for i in range(NoN):
#     plt.scatter(NL[i, 0], NL[i, 1], color='black')  # Вузли
#     plt.text(NL[i, 0], NL[i, 1], str(i + 1), color='red', fontsize=8)
#
#     # Побудова елементів
# for j in range(NoE):
#         nodes = EL[j]
#         x = NL[nodes, 0]
#         y = NL[nodes, 1]
#         x = np.append(x, x[0])  # Замкнути контур
#         y = np.append(y, y[0])
#         plt.plot(x, y, color='blue')
#         # Анотація елемента
#         cx = np.mean(x[:-1])
#         cy = np.mean(y[:-1])
#         plt.text(cx, cy, str(j + 1), color='green', fontsize=10)
#
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title(f'Сітка: {element_type}')
# plt.grid()
# plt.axis('equal')
# plt.show()
