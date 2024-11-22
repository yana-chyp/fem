import numpy as np

import base_functions_2d as bs2d
import finite_element_2d as fe2d
import system_2d as s2d
import nastia2d as n2d

def f(x, y):
    return 0

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
p = 4
m = 3

nodes, elements = n2d.uniform_mesh(d1, d2, p, m, 'D2QU4N')

# bounds = s2d.get_boundary_points(p, m)
# print(bounds)



matrix = s2d.set_up_matrix(d1, d2, p, m)
print(50*"-")
for row in matrix:
    print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )

print(50*"-")
f_vec = s2d.set_up_vector(f, base, d1, d2, p, m, element_type='D2QU4N')
print('[' + ', '.join([f"{el:.4f}" for el in f_vec]) + ']')

ug = [ug_3, ug_2, ug_4, ug_1]
matrix, f_vec = s2d.apply_boundary_conditions(matrix, f_vec, p, m, nodes, ug)

print(50*"-")
print(50*"-")

for row in matrix:
    print('[' + ', '.join([f"{el:.4f}" for el in row]) + ']' )
print(50*"-")

print('[' + ', '.join([f"{el:.4f}" for el in f_vec]) + ']')


res = np.linalg.solve(matrix, f_vec)
print('[' + ', '.join([f"{el:.4f}" for el in res]) + ']')

