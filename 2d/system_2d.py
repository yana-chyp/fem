# from distutils.command.install import value
import numpy as np

import mesh_2d as m2d
import base_functions_2d as bf2d
import finite_element_2d as fe2d
import base_functions_2d as bs2d
from base_functions_2d import ksi_right, ksi_left, eta_right, eta_left
import sympy as sp
import scipy.integrate as scin

def solve(b1, d1, b2, d2, p, m, degree, f, ug, element_type='D2QU4N'):
    nodes, elements = m2d.uniform_mesh(d1, d2, p, m, element_type, degree, b1, b2)
    h_x = (d1-b1)/p /degree; h_y = (d2-b2)/m /degree; J = h_x*h_y/((ksi_right-ksi_left)*(eta_right-eta_left))
    matrix = set_up_matrix(nodes, elements, degree, J, p, m)
    base = bs2d.get_base_functions(degree)
    f_vec = set_up_vector(f, base, nodes, elements, degree, b1, d1, b2, d2, p, m)
    matrix, f_vec = apply_boundary_conditions(matrix, f_vec, p, m, nodes, ug, degree)
    u = np.linalg.solve(matrix, f_vec)
    print('[' + ', '.join([f"{el:.4f}" for el in u]) + ']')
    m2d.plot_2d_solution(u, nodes, elements)


def set_up_matrix(nodes, elements, degree, J, p, m):
    # nodes, elements = m2d.uniform_mesh(d1, d2, p, m, element_type, degree)
    #npe - nodes per element (assume all are the same)
    npe = len(elements[0])
    #let it be rectangle

    #assume nodes are equidistant
    # h_x = (d1-b1)/p /degree; h_y = (d2-b2)/m /degree; J = h_x*h_y/((ksi_right-ksi_left)*(eta_right-eta_left))

    n = (degree*p+1)*(degree*m+1)
    matrix = [[0 for j in range(n)] for i in range(n)]

    em = fe2d.element_matrix(degree)
    #ioe - index of element
    for ioe in range(len(elements)):
        #noe - nodes of element
        noe = elements[ioe]
        # print(noe)
        for i in range(len(noe)):
            for j in range(len(noe)):
                im = noe[i]
                jm = noe[j]
                value = em[i][j]*J
                # print("(", i, ",  ", j, ") -> (", im, ", ", jm, "); value = ", f"{value:.4f}")
                matrix[im][jm] += value
            # print('[' + ', '.join([f"{el:.8f}" for el in matrix[i]]) + ']')
    return matrix


def set_up_vector(f, base, nodes, elements, degree, b1, d1, b2, d2, p, m):
    # nodes, elements = m2d.uniform_mesh(d1, d2, p, m, element_type, degree)

    n = (degree*p+1)*(m*degree+1)
    h_x = (d1-b1) / p /degree; h_y = (d2-b2) / m /degree
    # J = h_x * h_y / ((ksi_right - ksi_left) * (eta_right - eta_left))

    f_vec = [0 for i in range(n)]
    ksi = sp.symbols('ksi'); eta = sp.symbols('eta')
    x = sp.symbols('x'); y = sp.symbols('y')

    for ioe in range(len(elements)):
        noe = elements[ioe]
        x_0 = nodes[noe[0]][0]; y_0 = nodes[noe[0]][1]
        lin_ksi = 2*(x - x_0)/h_x - 1
        lin_eta = 2*(y - y_0)/h_y - 1
        f_expr = sp.sympify(f(x, y))
        for i in range(len(noe)):
            N_i = base[i].subs('ksi', lin_ksi)
            N_i = N_i.subs('eta', lin_eta)
            prod = sp.lambdify([x, y], N_i*f_expr)
            value = scin.dblquad(prod, x_0, x_0+h_x, y_0, y_0+h_y)[0]
            f_vec[noe[i]] += value
    return f_vec




def get_boundary_points(p, m, degree = 1):
    bounds = []
    gamma_1 = []; gamma_2 = []; gamma_3 = []; gamma_4 = []
    for j in range(degree*p+1):
        gamma_1.append(j)
    bounds.append(gamma_1)
    for i in range(1, degree*m+1):
        index = (i+1)*(degree*p+1)-1
        gamma_2.append(index)
    bounds.append(gamma_2)
    for j in reversed(range(degree*p)):
        index = (degree*p+1)*degree*m + j
        gamma_3.append(index)
    bounds.append(gamma_3)
    for i in reversed(range(1, degree*m)):
        index = i*(degree*p+1)
        gamma_4.append(index)
    bounds.append(gamma_4)
    return bounds

def apply_boundary_conditions(matrix, f_vec, p, m, nodes, ug, degree = 1):
    bounds = get_boundary_points(p, m, degree)
    # print(bounds)
    n = len(ug)
    for i in range(n):
        for index in bounds[i]:
            for j in range(len(matrix[index])):
                matrix[index][j] = 0
            matrix[index][index] = 1
            f_vec[index] = ug[i](nodes[index][0], nodes[index][1])
    return (matrix, f_vec)






