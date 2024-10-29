import base_functions as base
import finite_element as finel
import nastia as nastia
from fem_lib import get_base_function
import matplotlib.pyplot as plt
import sympy as sp

def f(ksi):
    # if i==1:
    return (ksi+1/3)*(ksi-1/3)*(ksi-1)/(- 16/9)

    # return 0

a = -1
b = 1
n = 5
m = 3
mesh = nastia.get_mesh(a, b, n)

ksi_element = finel.element_matrix(m)
for i in range(m + 1):
    print(ksi_element[i])

# div = 100
#
# h = (base.ksi_right - base.ksi_left) / div
# ksi_s = [base.ksi_left + k*h for k in range(div+1)]
# # print(ksi_s)
# for i in range(n) :
#     base_functions = base.get_base_functions(m)
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


