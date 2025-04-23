import numpy as np
import sympy as sp
import mesh_2d as m2d
import system_2d as s2d
from enum import Enum
import matplotlib.pyplot as plt


# def f(x, y):
#      return 10 if ((x - 0.5)**2 + (y - 0.5)**2) < 0.04 else 0


#     return x* (1-x)* y *(1-y)
#
#
#
# def exact_solution(x, y):
#     return x* (1-x)* y *(1-y)

# def ug_1(x, y):
#     return 0
# def ug_2(x, y):
#     return (-(y-0.5)**2 + 0.25)*8 + 2
# def ug_3(x, y):
#     return x
# def ug_4(x, y):
#     return x
# def ug_1(x, y):  # верх
#     return 1
# def ug_2(x, y):  # низ
#     return 0
# def ug_3(x, y):  # ліва
#     return 0  # Нейман
# def ug_4(x, y):  # права
#     return 0  # Нейман
#
# def k(x, y):
#     return 1 + x**2 + y**2

#     return 2;


def N(i, ksi, eta, ap):
    if ap == 1:
        functions = [
            (1 - ksi) * (1 - eta) / 4,
            (1 + ksi) * (1 - eta) / 4,
            (1 + ksi) * (1 + eta) / 4,
            (1 - ksi) * (1 + eta) / 4
        ]
    elif ap == 2:
        functions = [
            (1 / 4) * (1 - ksi) * (1 - eta) * ksi * eta,
            (1 / 4) * (ksi + 1) * (1 - eta) * ksi * ( - eta),
            (1 / 4) * (ksi + 1) * (eta + 1) * ksi * eta,
            (1 / 4) * (1 - ksi) * (eta + 1) * ( - ksi) * eta,
            (1 / 2) * (1 - ksi ** 2) * ( - eta ) * (1 - eta),
            (1 / 2) * (1 - eta ** 2) * ksi * (ksi + 1),
            (1 / 2) * (1 - ksi ** 2) * eta * (eta + 1),
            (1 / 2) * (1 - eta ** 2) * ( - ksi ) * (1 - ksi),
            (1 - ksi ** 2) * (1 - eta ** 2)
        ]
    elif ap == 3:
        functions = [
            (81 / 256) * (1 - ksi) * (1 - eta) * ((1 / 9) - ksi ** 2) * ((1 / 9) - eta ** 2),
            (81 / 256) * (1 + ksi) * (1 - eta) * ((1 / 9) - ksi ** 2) * ((1 / 9) - eta ** 2),
            (81 / 256) * (1 + ksi) * (1 + eta) * ((1 / 9) - ksi ** 2) * ((1 / 9) - eta ** 2),
            (81 / 256) * (1 - ksi) * (1 + eta) * ((1 / 9) - ksi ** 2) * ((1 / 9) - eta ** 2),

            (243 / 256) * (1 - ksi ** 2) * (eta ** 2 - (1 / 9)) * ((1 / 3) - ksi) * (1 - eta),
            (243 / 256) * (1 - eta ** 2) * (ksi ** 2 - (1 / 9)) * ((1 / 3) - eta) * (1 + ksi),
            (243 / 256) * (1 - ksi ** 2) * (eta ** 2 - (1 / 9)) * ((1 / 3) + ksi) * (1 + eta),
            (243 / 256) * (1 - eta ** 2) * (ksi ** 2 - (1 / 9)) * ((1 / 3) + eta) * (1 - ksi),
            (243 / 256) * (1 - ksi ** 2) * (eta ** 2 - (1 / 9)) * ((1 / 3) + ksi) * (1 - eta),
            (243 / 256) * (1 - eta ** 2) * (ksi ** 2 - (1 / 9)) * ((1 / 3) + eta) * (1 + ksi),
            (243 / 256) * (1 - ksi ** 2) * (eta ** 2 - (1 / 9)) * ((1 / 3) - ksi) * (1 + eta),
            (243 / 256) * (1 - eta ** 2) * (ksi ** 2 - (1 / 9)) * ((1 / 3) - eta) * (1 - ksi),

            (81 / 256) * (1 - ksi ** 2) * (1 - 3 * ksi) * (1 - eta ** 2) * (1 - 3 * eta),
            (81 / 256) * (1 - ksi ** 2) * (1 + 3 * ksi) * (1 - eta ** 2) * (1 - 3 * eta),
            (81 / 256) * (1 - ksi ** 2) * (1 + 3 * ksi) * (1 - eta ** 2) * (1 + 3 * eta),
            (81 / 256) * (1 - ksi ** 2) * (1 - 3 * ksi) * (1 - eta ** 2) * (1 + 3 * eta)
        ]
    return functions[i]

def isoparametric_transform(ksi, eta, x_coords, y_coords, ap):
    num_nodes = 4 if ap == 1 else (9 if ap == 2 else 16)
    x = sum(N(i, ksi, eta, ap) * x_coords[i] for i in range(num_nodes))
    y = sum(N(i, ksi, eta, ap) * y_coords[i] for i in range(num_nodes))
    return x, y



def uniform_mesh_with_vertices(vertices, p, m, element_type, degree=1):
    x0, y0 = vertices[0]
    x1, y1 = vertices[1]
    x2, y2 = vertices[2]
    x3, y3 = vertices[3]

    nodes_x = degree * p + 1
    nodes_y = degree * m + 1

    NL = np.zeros([nodes_x * nodes_y, 2])

    if degree == 1:
        x_coords = np.array([x0, x1, x2, x3])
        y_coords = np.array([y0, y1, y2, y3])
    elif degree == 2:
        x4, y4 = (x1 + x0) / 2, (y1 + y0) / 2
        x5, y5 = (x1 + x2) / 2, (y1 + y2) / 2
        x6, y6 = (x3 + x2) / 2, (y3 + y2) / 2
        x7, y7 = (x3 + x0) / 2, (y3 + y0) / 2
        x8, y8 = (x1 + x2 + x3 + x0) / 4, (y1 + y2 + y3 + y0) / 4
        x_coords = np.array([x0, x1, x2, x3, x4, x5, x6, x7, x8])
        y_coords = np.array([y0, y1, y2, y3, y4, y5, y6, y7, y8])
    elif degree == 3:
        x4, y4 = (2 * x0 + x1) / 3, (2 * y0 + y1) / 3
        x5, y5 = (x0 + 2 * x1) / 3, (y0 + 2 * y1) / 3
        x6, y6 = (2 * x1 + x2) / 3, (2 * y1 + y2) / 3
        x7, y7 = (x1 + 2 * x2) / 3, (y1 + 2 * y2) / 3
        x8, y8 = (2 * x2 + x3) / 3, (2 * y2 + y3) / 3
        x9, y9 = (x2 + 2 * x3) / 3, (y2 + 2 * y3) / 3
        x10, y10 = (2 * x3 + x0) / 3, (2 * y3 + y0) / 3
        x11, y11 = (x3 + 2 * x0) / 3, (y3 + 2 * y0) / 3
        x12, y12 = (2 * x11 + x6) / 3, (2 * y11 + y6) / 3
        x13, y13 = (x11 + 2 * x6) / 3, (y11 + 2 * y6) / 3
        x14, y14 = (x10 + 2 * x7) / 3, (y10 + 2 * y7) / 3
        x15, y15 = (2 * x10 + x7) / 3, (2 * y10 + y7) / 3
        x_coords = np.array([x0, x1, x2, x3, x4, x6, x8, x10, x5, x7, x9, x11, x12, x13, x14, x15])
        y_coords = np.array([y0, y1, y2, y3, y4, y6, y8, y10, y5, y7, y9, y11, y12, y13, y14, y15])

    for i in range(nodes_y):
        for j in range(nodes_x):
            ksi = 2 * j / (nodes_x - 1) - 1
            eta = 2 * i / (nodes_y - 1) - 1

            x, y = isoparametric_transform(ksi, eta, x_coords, y_coords, degree)
            NL[i * nodes_x + j, :] = [x, y]

    if degree == 1:
        NPE = 4
    elif degree == 2:
        NPE = 9
    elif degree == 3:
        NPE = 16
    else:
        raise ValueError("Параметр degree повинен бути 1, 2 або 3.")

    NoE = p * m
    EL = np.zeros([NoE, NPE], dtype=int)

    e = 0
    for i in range(m):
        for j in range(p):
            base = i * degree * nodes_x + j * degree
            if degree == 1:
                # 4 вузли на елемент
                n1 = base
                n2 = n1 + 1
                n3 = n1 + nodes_x
                n4 = n3 + 1
                EL[e, :] = [n1, n2, n4, n3]
            elif degree == 2:
                # 9 вузлів на елемент
                n1 = base #1
                n2 = n1 + 2 #2
                n3 = n1 + 2 * nodes_x #4
                n4 = n3 + 2 #3
                n5 = n1 + 1 #5
                n6 = n3 + 1 #7
                n7 = n1 + nodes_x #8
                n8 = n7 + 2 #6
                n9 = n7 + 1 #9
                EL[e, :] = [n1, n2, n4, n3, n5, n8, n6, n7, n9]
            elif degree == 3:
                # 16 вузлів на елемент
                n1 = base #1
                n2 = n1 + 3  #2
                n3 = n1 + 3 * nodes_x #4
                n4 = n3 + 3 #3
                n5 = n1 + 1 #5
                n6 = n1 + 2 #9
                n7 = n3 + 1 #11
                n8 = n3 + 2 #7
                n9 = n1 + nodes_x #12
                n10 = n9 + 1 #13
                n11 = n9 + 2 #14
                n12 = n9 + 3 #6
                n13 = n1 + 2 * nodes_x #8
                n14 = n13 + 1 #16
                n15 = n13 + 2 #15
                n16 = n13 + 3 #10
                EL[e, :] = [n1, n2, n4, n3, n5, n12, n8, n13, n6, n16, n7, n9, n10, n11, n15, n14]

            e += 1

    return NL, EL



def set_up_vector(f, nodes, elements,p, m, ap):
    f_vec = np.zeros((ap*p+1)*(ap*m+1))

    dN_dksi_list, dN_deta_list = compute_partial_derivatives(ap)

    if ap == 1:
        gauss_points = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
        gauss_weights = [1.0, 1.0]

    elif ap == 2:
        gauss_points = [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5)]
        gauss_weights = [5 / 9, 8 / 9, 5 / 9]

    elif ap == 3:
        gauss_points = [
            -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
            -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
            np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
            np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5))
        ]
        gauss_weights = [
            (18 - np.sqrt(30)) / 36,
            (18 + np.sqrt(30)) / 36,
            (18 - np.sqrt(30)) / 36,
            (18 + np.sqrt(30)) / 36
        ]

    # Обчислення інтегралів для кожного елемента
    for ioe, noe in enumerate(elements):
        x_coords = [nodes[i][0] for i in noe]
        y_coords = [nodes[i][1] for i in noe]

        for i in range(len(noe)):
            integral_value = 0.0

            for ksi_idx, ksi_point in enumerate(gauss_points):
                for eta_idx, eta_point in enumerate(gauss_points):
                    x_val, y_val = isoparametric_transform(ksi_point, eta_point, x_coords, y_coords, ap)
                    J_val = abs(compute_jacobian(ksi_point, eta_point, x_coords, y_coords, dN_dksi_list, dN_deta_list))
                    f_val = f(x_val, y_val)
                    N_i = N(i, ksi_point, eta_point, ap)

                    # Використання відповідних ваг Гаусса
                    integral_value += N_i * f_val * J_val * gauss_weights[ksi_idx] * gauss_weights[eta_idx]

            f_vec[noe[i]] += integral_value

    return f_vec

def compute_element_stiffness(elements, nodes, ap, k1, k2):
    ksi, eta = sp.symbols('ksi eta')
    element_matrices = []

    dN_dksi_list, dN_deta_list = compute_partial_derivatives(ap)

    if ap == 1:
        gauss_points = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
        gauss_weights = [1.0, 1.0]
    elif ap == 2:
        gauss_points = [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5)]
        gauss_weights = [5 / 9, 8 / 9, 5 / 9]
    elif ap == 3:
        gauss_points = [
            -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
            -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
            np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
            np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5))
        ]
        gauss_weights = [
            (18 - np.sqrt(30)) / 36,
            (18 + np.sqrt(30)) / 36,
            (18 - np.sqrt(30)) / 36,
            (18 + np.sqrt(30)) / 36
        ]

    for noe in elements:
        x_coords = [nodes[i][0] for i in noe]
        y_coords = [nodes[i][1] for i in noe]

        K_local = np.zeros((len(noe), len(noe)))

        for i in range(len(noe)):
            for j in range(len(noe)):
                integral_value = 0.0

                for ksi_idx, ksi_point in enumerate(gauss_points):
                    for eta_idx, eta_point in enumerate(gauss_points):
                        J_val = compute_jacobian(ksi_point, eta_point, x_coords, y_coords, dN_dksi_list, dN_deta_list)
                        J_inv = 1 / J_val


                        dN_i_dksi = dN_dksi_list[i](ksi_point, eta_point)
                        dN_i_deta = dN_deta_list[i](ksi_point, eta_point)
                        dN_j_dksi = dN_dksi_list[j](ksi_point, eta_point)
                        dN_j_deta = dN_deta_list[j](ksi_point, eta_point)


                        x_gp = sum(x_coords[n] * dN_dksi_list[n](ksi_point, eta_point) for n in range(len(noe)))
                        y_gp = sum(y_coords[n] * dN_deta_list[n](ksi_point, eta_point) for n in range(len(noe)))

                        # Будування матриці
                        # Значення коефіцієнтів
                        k1_val = k1(x_gp, y_gp)
                        k2_val = k2(x_gp, y_gp)

                        # Інтегральний доданок з коефіцієнтом k(x, y)
                        integrand = (k1_val * dN_i_dksi * dN_j_dksi + k2_val * dN_i_deta * dN_j_deta) * J_inv
                        integral_value += integrand * gauss_weights[ksi_idx] * gauss_weights[eta_idx]

                K_local[i, j] = integral_value

        element_matrices.append(K_local)

    return element_matrices

# def compute_element_stiffness(elements, nodes, ap, k):
#     ksi, eta = sp.symbols('ksi eta')
#     element_matrices = []
#
#     dN_dksi_list, dN_deta_list = compute_partial_derivatives(ap)
#
#     if ap == 1:
#         gauss_points = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
#         gauss_weights = [1.0, 1.0]
#     elif ap == 2:
#         gauss_points = [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5)]
#         gauss_weights = [5 / 9, 8 / 9, 5 / 9]
#     elif ap == 3:
#         gauss_points = [
#             -np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
#             -np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5)),
#             np.sqrt((3 / 7) - (2 / 7) * np.sqrt(6 / 5)),
#             np.sqrt((3 / 7) + (2 / 7) * np.sqrt(6 / 5))
#         ]
#         gauss_weights = [
#             (18 - np.sqrt(30)) / 36,
#             (18 + np.sqrt(30)) / 36,
#             (18 - np.sqrt(30)) / 36,
#             (18 + np.sqrt(30)) / 36
#         ]
#
#     for noe in elements:
#         x_coords = [nodes[i][0] for i in noe]
#         y_coords = [nodes[i][1] for i in noe]
#
#         K_local = np.zeros((len(noe), len(noe)))
#
#         for i in range(len(noe)):
#             for j in range(len(noe)):
#                 integral_value = 0.0
#
#                 for ksi_idx, ksi_point in enumerate(gauss_points):
#                     for eta_idx, eta_point in enumerate(gauss_points):
#                         J_val = compute_jacobian(ksi_point, eta_point, x_coords, y_coords, dN_dksi_list, dN_deta_list)
#                         J_inv = 1 / J_val
#
#                         # dN_i_dksi = dN_dksi_list[i].subs({ksi: ksi_point, eta: eta_point}).evalf()
#                         # dN_i_deta = dN_deta_list[i].subs({ksi: ksi_point, eta: eta_point}).evalf()
#                         # dN_j_dksi = dN_dksi_list[j].subs({ksi: ksi_point, eta: eta_point}).evalf()
#                         # dN_j_deta = dN_deta_list[j].subs({ksi: ksi_point, eta: eta_point}).evalf()
#                         dN_i_dksi = dN_dksi_list[i](ksi_point, eta_point)
#                         dN_i_deta = dN_deta_list[i](ksi_point, eta_point)
#                         dN_j_dksi = dN_dksi_list[j](ksi_point, eta_point)
#                         dN_j_deta = dN_deta_list[j](ksi_point, eta_point)
#
#                         # integrand = (dN_i_dksi * dN_j_dksi + dN_i_deta * dN_j_deta) * J_inv
#                         # integral_value += integrand * gauss_weights[ksi_idx] * gauss_weights[eta_idx]
#                         # Обчислюємо фізичні координати точки інтегрування
#                         # x_gp = sum(x_coords[n] * float(dN_dksi_list[n].subs({ksi: ksi_point, eta: eta_point})) for n in
#                         #            range(len(noe)))
#                         # y_gp = sum(y_coords[n] * float(dN_deta_list[n].subs({ksi: ksi_point, eta: eta_point})) for n in
#                         #            range(len(noe)))
#                         x_gp = sum(x_coords[n] * dN_dksi_list[n](ksi_point, eta_point) for n in range(len(noe)))
#                         y_gp = sum(y_coords[n] * dN_deta_list[n](ksi_point, eta_point) for n in range(len(noe)))
#
#                         # Будування матриці
#                         # Значення коефіцієнта k в цій точці
#                         k_val = k(x_gp, y_gp)
#
#                         # Інтегральний доданок з коефіцієнтом k(x, y)
#                         integrand = k_val * (dN_i_dksi * dN_j_dksi + dN_i_deta * dN_j_deta) * J_inv
#                         integral_value += integrand * gauss_weights[ksi_idx] * gauss_weights[eta_idx]
#
#                 K_local[i, j] = integral_value
#
#         element_matrices.append(K_local)
#
#     return element_matrices


def assemble_global_stiffness_matrix(nodes, elements, p, m, element_matrices, ap):
    num_nodes = (ap*p+1)*(ap*m+1)
    K = np.zeros((num_nodes, num_nodes))

    for elem_idx, noe in enumerate(elements):
        em = element_matrices[elem_idx]

        for i in range(len(noe)):
            for j in range(len(noe)):
                im = noe[i]
                jm = noe[j]
                K[im, jm] += em[i, j]

    return K

# def compute_jacobian(ksi, eta, x_coords, y_coords, dN_dksi_list, dN_deta_list):
#     dx_dksi = sum(dN_dksi_list[j].subs({sp.symbols('ksi'): ksi, sp.symbols('eta'): eta}).evalf() * x_coords[j] for j in range(len(x_coords)))
#     dx_deta = sum(dN_deta_list[j].subs({sp.symbols('ksi'): ksi, sp.symbols('eta'): eta}).evalf() * x_coords[j] for j in range(len(x_coords)))
#     dy_dksi = sum(dN_dksi_list[j].subs({sp.symbols('ksi'): ksi, sp.symbols('eta'): eta}).evalf() * y_coords[j] for j in range(len(y_coords)))
#     dy_deta = sum(dN_deta_list[j].subs({sp.symbols('ksi'): ksi, sp.symbols('eta'): eta}).evalf() * y_coords[j] for j in range(len(y_coords)))
#
#     J_val = abs(dx_dksi * dy_deta - dx_deta * dy_dksi)
#     return J_val
def compute_jacobian(ksi, eta, x_coords, y_coords, dN_dksi_funcs, dN_deta_funcs):
    dx_dksi = sum(dN_dksi_funcs[j](ksi, eta) * x_coords[j] for j in range(len(x_coords)))
    dx_deta = sum(dN_deta_funcs[j](ksi, eta) * x_coords[j] for j in range(len(x_coords)))
    dy_dksi = sum(dN_dksi_funcs[j](ksi, eta) * y_coords[j] for j in range(len(y_coords)))
    dy_deta = sum(dN_deta_funcs[j](ksi, eta) * y_coords[j] for j in range(len(y_coords)))

    J_val = abs(dx_dksi * dy_deta - dx_deta * dy_dksi)
    return J_val

# Функція для обчислення частинних похідних
# def compute_partial_derivatives(ap):
#     ksi, eta = sp.symbols('ksi eta')
#
#     num_nodes = 4 if ap == 1 else (9 if ap == 2 else 16)
#     dN_dksi_list = []
#     dN_deta_list = []
#
#     for i in range(num_nodes):
#         Ni = N(i, ksi, eta, ap)
#         dN_dksi = sp.diff(Ni, ksi)
#         dN_deta = sp.diff(Ni, eta)
#         dN_dksi_list.append(dN_dksi)
#         dN_deta_list.append(dN_deta)
#
#     return dN_dksi_list, dN_deta_list
def compute_partial_derivatives(ap):
    ksi, eta = sp.symbols('ksi eta')

    num_nodes = 4 if ap == 1 else (9 if ap == 2 else 16)
    dN_dksi_funcs = []
    dN_deta_funcs = []

    for i in range(num_nodes):
        Ni = N(i, ksi, eta, ap)
        dN_dksi = sp.diff(Ni, ksi)
        dN_deta = sp.diff(Ni, eta)

        # Використовуємо lambdify замість подальших підстановок у циклі
        dN_dksi_func = sp.lambdify((ksi, eta), dN_dksi, modules="numpy")
        dN_deta_func = sp.lambdify((ksi, eta), dN_deta, modules="numpy")

        dN_dksi_funcs.append(dN_dksi_func)
        dN_deta_funcs.append(dN_deta_func)

    return dN_dksi_funcs, dN_deta_funcs


def apply_boundary_conditions(matrix, f_load, p, m, nodes, ug, ap = 1):
    validate_boundary_conditions(ug)
    bounds = get_boundary_elements_and_nodes(p, m, ug, ap)
    n = len(ug)
    for i in range(n):
        for index in bounds[i][1]:  # [1] — список вузлів
            for j in range(len(matrix[index])):
                matrix[index][j] = 0
            matrix[index][index] = 1
            f_load[index] = ug[i][1](nodes[index][0], nodes[index][1])  # <== ось тут виправлено
    return (matrix, f_load)


def get_boundary_elements_and_nodes(p, m, ug, degree = 1):
    bounds = []
    #gamma_i consists of [elements], [nodes]
    elements_1 = [j for j in range(p)]
    nodes_1 = [j for j in range(1, degree*p)]

    elements_2 = [p*(i+1)-1 for i in range(m)]
    nodes_2 = [(i+1)*(degree*p+1)-1 for i in range(1, degree*m)]

    intersection = degree*p
    # print(intersection)
    if (ug[0][0]==TypeOfBoundCond.DIRICHLET):
        nodes_1.append(intersection)
        # print('1', nodes_1)
    else:
        nodes_2.insert(0, intersection)
        # print('2', nodes_2)

    elements_3 = [p*m-1-j for j in range(p)]
    nodes_3 = [(degree*p+1)*degree*m + j for j in reversed(range(1, degree*p))]

    intersection = (degree*p+1)*(degree*m+1)-1
    # print(intersection)
    if (ug[1][0]==TypeOfBoundCond.DIRICHLET):
        nodes_2.append(intersection)
        # print('2', nodes_2)
    else:
        nodes_3.insert(0, intersection)
        # print('3', nodes_3)

    elements_4 = [p*i for i in reversed(range(m))]
    nodes_4 = [i*(degree*p+1) for i in reversed(range(1, degree*m))]

    intersection = degree*m*(degree*p+1)
    # print(intersection)
    if (ug[2][0]==TypeOfBoundCond.DIRICHLET):
        nodes_3.append(intersection)
        # print('3', nodes_3)
    else:
        nodes_4.insert(0, intersection)
        # print('4', nodes_4)

    intersection = 0
    # print(intersection)
    if (ug[3][0]==TypeOfBoundCond.DIRICHLET):
        nodes_4.append(intersection)
        # print('4', nodes_4)
    else:
        nodes_1.insert(0, intersection)
        # print('1', nodes_1)

    gamma_1 = [elements_1, nodes_1]
    bounds.append(gamma_1)
    gamma_2 = [elements_2, nodes_2]
    bounds.append(gamma_2)
    gamma_3 = [elements_3, nodes_3]
    bounds.append(gamma_3)
    gamma_4 = [elements_4, nodes_4]
    bounds.append(gamma_4)
    return bounds


class TypeOfBoundCond(Enum):
    DIRICHLET = 1
    NEUMANN = 2


def validate_boundary_conditions(ug):
    isPresentDirichlet = False
    for ug_ in ug:
        if ug_[0] == TypeOfBoundCond.DIRICHLET:
            isPresentDirichlet = True
            break
    if not isPresentDirichlet:
        raise Exception("Wrong boundary conditions: all are Neumann, Dirichlet must be present")

def get_boundary_points(p, m, ap = 1):
    bounds = []
    gamma_1 = []; gamma_2 = []; gamma_3 = []; gamma_4 = []
    for j in range(ap*p+1):
        gamma_1.append(j)
    bounds.append(gamma_1)
    for i in range(1, ap*m+1):
        index = (i+1)*(ap*p+1)-1
        gamma_2.append(index)
    bounds.append(gamma_2)
    for j in reversed(range(ap*p)):
        index = (ap*p+1)*ap*m + j
        gamma_3.append(index)
    bounds.append(gamma_3)
    for i in reversed(range(1, ap*m)):
        index = i*(ap*p+1)
        gamma_4.append(index)
    bounds.append(gamma_4)
    return bounds

# # Основна функція
# def main():
#     #vertices = [(0, 0), (1, 0.25), (1, 0.75), (0, 1)]
#     vertices = [(0, 0), (1, 0), (1, 1), (0, 1)]
#     element_type = 'D2QU4N'
#     p = 5  # Кількість елементів по x
#     m = 5  # Кількість елементів по y
#     ap = 1
#     # ug = [ug_3, ug_2, ug_4, ug_1]
#     ug = [[TypeOfBoundCond.NEUMANN, ug_3],
#           [TypeOfBoundCond.DIRICHLET, ug_1],
#           [TypeOfBoundCond.NEUMANN, ug_4],
#           [TypeOfBoundCond.DIRICHLET, ug_2]]
#
#     NL, EL= uniform_mesh_with_vertices(vertices, p, m, element_type, ap)
#     f_load = set_up_vector(f, NL, EL, p, m, ap)
#     elem_matrices = compute_element_stiffness(EL, NL, ap, k)
#     matrix =  assemble_global_stiffness_matrix(NL, EL, p, m, elem_matrices, ap)
#     print(matrix)
#     matrix, f_load = apply_boundary_conditions(matrix, f_load, p, m, NL, ug, ap)
#     u = np.linalg.solve(matrix, f_load)
#
#     m2d.plot_2d_solution(u, NL, EL)
#     print("Координати вузлів:\n", NL)
#     print("Елементи:\n", EL)
#     print(f_load)
#     print(matrix)
#     print(u)
#
#     errors = []
#     hs = []
#
#     for N in [4, 8, 16, 32]:
#         p = m = N
#         NL, EL = uniform_mesh_with_vertices(vertices, p, m, element_type, ap)
#         f_load = set_up_vector(f, NL, EL, p, m, ap)
#         elem_matrices = compute_element_stiffness(EL, NL, ap, k)
#         matrix = assemble_global_stiffness_matrix(NL, EL, p, m, elem_matrices, ap)
#         matrix, f_load = apply_boundary_conditions(matrix, f_load, p, m, NL, ug, ap)
#         u = np.linalg.solve(matrix, f_load)
#
#         # Створення базисних функцій
#         base_functions = [lambda ξ, η, i=i: N(i, ξ, η, ap) for i in range(len(EL[0]))]
#
#         # Знаходження розміру елемента (h)
#         h = 1.0 / N
#         hs.append(h)
#
#         # Обчислення похибки
#         error = calculate_L2_error(u_exact, u, NL, EL, base_functions, degree=ap + 1,
#                                    J=0.25)  # J може залежати від елемента!
#         errors.append(error)

def plot_statistics(start, end, statistics, u_exact, n_points=100):
    """
    Малює u вздовж відрізка від точки start до end.
    """
    abscissas = np.linspace(0, 1, n_points)
    xs = [(1 - t) * start[0] + t * end[0] for t in abscissas]
    ys = [(1 - t) * start[1] + t * end[1] for t in abscissas]
    u_values = [u_exact(x, y) for x, y in zip(xs, ys)]

    plt.plot(abscissas, u_values, color='red', label='u exact')

    for i, (x_approx, u_approx) in enumerate(statistics):
        plt.plot(x_approx, u_approx, label=f'approx {i}')

    plt.xlabel('normalized distance along line')
    plt.ylabel('u')
    plt.title('Exact vs Approx Solutions')
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_errors(errors):
    plt.plot([i for i in range(len(errors))], errors, color='green')
    plt.xlabel('degree of nodes')
    plt.ylabel('errors')
    plt.grid(True)
    plt.show()


def interpolate_solution(x, y, u, NL, EL, ap):
    from scipy.optimize import root

    num_nodes = 4 if ap == 1 else (9 if ap == 2 else 16)

    for element in EL:
        nodes_coords = [NL[i] for i in element]
        x_coords = [pt[0] for pt in nodes_coords]
        y_coords = [pt[1] for pt in nodes_coords]

        def equations(p):
            ksi, eta = p
            x_mapped = sum(N(i, ksi, eta, ap) * x_coords[i] for i in range(num_nodes))
            y_mapped = sum(N(i, ksi, eta, ap) * y_coords[i] for i in range(num_nodes))
            return [x_mapped - x, y_mapped - y]

        sol = root(equations, [0, 0])

        if sol.success:
            ksi, eta = sol.x
            if -1 <= ksi <= 1 and -1 <= eta <= 1:
                u_local = [u[i] for i in element]
                u_val = sum(u_local[i] * N(i, ksi, eta, ap) for i in range(num_nodes))
                return u_val

    return 0  # якщо точка не належить жодному елементу


# def u_exact(x, y):
#     return np.sin(np.pi * x) * np.sin(np.pi * y)
#
#
# def f(x, y):
#     return 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
#
#
# def ug_1(x, y):  # верхнє ребро (y = 1)
#     return np.sin(np.pi * x) * np.sin(np.pi * 1)
#
# def ug_2(x, y):  # ліве ребро (x = 0)
#     return np.sin(np.pi * 0) * np.sin(np.pi * y)
#
# def ug_3(x, y):  # нижнє ребро (y = 0)
#     return np.sin(np.pi * x) * np.sin(np.pi * 0)
#
# def ug_4(x, y):  # праве ребро (x = 1)
#     return np.sin(np.pi * 1) * np.sin(np.pi * y)

def u_exact(x, y):
    return y**2
    # return x**2+y
    #return y**2 + x
    # return x**2 - y**2

def f(x, y):
    return 2

def ug_1(x, y):
    return y**2
    # return y
    # return y**2
def ug_2(x, y):
    return y**2
    # return y+1
    # return y**2 + 1
def ug_3(x, y):
    return 0
    # return x**2
    # return x
def ug_4(x, y):
    return 1
    # return x**2+1
    # return x+1

# def k(x, y):
#     return 1


def calculate_L2_error(u_exact, u_values, nodes, elements, base, degree):
    from numpy.polynomial.legendre import leggauss
    dN_dksi_list, dN_deta_list = compute_partial_derivatives(degree)
    nq = degree + 1  # кількість точок для квадратури
    quad_points, quad_weights = leggauss(nq)  # квадратура на відрізку [-1, 1]

    error_squared = 0

    for element in elements:
        # координати вузлів елемента
        element_coords = [nodes[i] for i in element]
        # значення u_h у вузлах цього елемента
        u_local = [u_values[i] for i in element]

        # проходимо по всім точкам квадратури (ξ, η)
        for i in range(nq):
            for j in range(nq):
                ξ = quad_points[i]
                η = quad_points[j]
                w = quad_weights[i] * quad_weights[j]

                # обчислюємо x, y через базис і координати вузлів
                x = sum(base[k](ξ, η) * element_coords[k][0] for k in range(len(element)))
                y = sum(base[k](ξ, η) * element_coords[k][1] for k in range(len(element)))

                u_h_val = sum(u_local[k] * base[k](ξ, η) for k in range(len(element)))
                u_ex_val = u_exact(x, y)

                J_val = compute_jacobian(ξ, η, [coord[0] for coord in element_coords],[coord[1] for coord in element_coords], dN_dksi_list, dN_deta_list)
                error_squared += (u_ex_val - u_h_val) ** 2 * w * J_val
    return np.sqrt(error_squared)

def main():
    vertices = [(0, 0), (1, 0), (1, 1), (0, 1)]
    k1 = lambda x, y: 2.0  # по x
    k2 = lambda x, y: 1.0  # по y
    ap = 1

    ug = [
        [TypeOfBoundCond.DIRICHLET, ug_3],  # нижнє (y = 0)
        [TypeOfBoundCond.DIRICHLET, ug_2],  # ліве (x = 0)
        [TypeOfBoundCond.DIRICHLET, ug_4],  # праве (x = 1)
        [TypeOfBoundCond.DIRICHLET, ug_1]  # верхнє (y = 1)
    ]

    # vertices = [(0, 0), (1, 0.25), (1, 0.75), (0, 1)]
    #     vertices = [(0, 0), (1, 0), (1, 1), (0, 1)]
    element_type = 'D2QU4N'
    p = 5  # Кількість елементів по x
    m = 5  # Кількість елементів по y
    # # ug = [ug_3, ug_2, ug_4, ug_1]
    # # ug = [[TypeOfBoundCond.NEUMANN, ug_3],
    # #       [TypeOfBoundCond.DIRICHLET, ug_1],
    # #       [TypeOfBoundCond.NEUMANN, ug_4],
    # #       [TypeOfBoundCond.DIRICHLET, ug_2]]
    #
    NL, EL= uniform_mesh_with_vertices(vertices, p, m, element_type, ap)
    f_load = set_up_vector(f, NL, EL, p, m, ap)
    elem_matrices = compute_element_stiffness(EL, NL, ap, k1=k1, k2=k2)
    matrix =  assemble_global_stiffness_matrix(NL, EL, p, m, elem_matrices, ap)
    print(matrix)
    matrix, f_load = apply_boundary_conditions(matrix, f_load, p, m, NL, ug, ap)
    u = np.linalg.solve(matrix, f_load)

    m2d.plot_2d_solution(u, NL, EL)
    print("Координати вузлів:\n", NL)
    print("Елементи:\n", EL)
    print(f_load)
    print(matrix)
    print(u)
#

    # errors = []
    # hs = []
    # statistics = []
    # start = (0, 0)
    # end = (1, 1)
    # # start = (0, 0.5)
    # # end = (1, 0.5)
    #
    # for num in [4, 8, 16, 32]:
    #     p = m = num
    #     NL, EL = uniform_mesh_with_vertices(vertices, p, m, element_type, ap)
    #     f_load = set_up_vector(f, NL, EL, p, m, ap)
    #     elem_matrices = compute_element_stiffness(EL, NL, ap, k)
    #     matrix = assemble_global_stiffness_matrix(NL, EL, p, m, elem_matrices, ap)
    #     matrix, f_load = apply_boundary_conditions(matrix, f_load, p, m, NL, ug, ap)
    #     u = np.linalg.solve(matrix, f_load)
    #     print("Максимальне u:", max(u))
    #     print("Мінімальне u:", min(u))
    #     x_line = []
    #     u_line = []
    #     for t in np.linspace(0, 1, 100):
    #         x = (1 - t) * start[0] + t * end[0]
    #         y = (1 - t) * start[1] + t * end[1]
    #         u_val = interpolate_solution(x, y, u, NL, EL, ap)
    #         x_line.append(t)
    #         u_line.append(u_val)
    #     statistics.append((x_line, u_line))
    #
    #     # Створення базисних функцій
    #     base_functions = [lambda ξ, η, i=i: N(i, ξ, η, ap) for i in range(len(EL[0]))]
    #
    #     # Знаходження розміру елемента (h)
    #     h = 1.0 / num
    #     hs.append(h)
    #
    #     # Обчислення похибки
    #     error = calculate_L2_error(u_exact, u, NL, EL, base_functions, degree=ap + 1)
    #     errors.append(error)
    #
    # # Візуалізація збіжності
    # plot_statistics(start, end, statistics, u_exact)
    #
    # plot_errors(errors)


if __name__ == "__main__":
    main()
