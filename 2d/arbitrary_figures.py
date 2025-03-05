import numpy as np
import sympy as sp
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
    return x
def ug_4(x, y):
    return x


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


def compute_element_stiffness(elements, nodes, ap):
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

                        dN_i_dksi = dN_dksi_list[i].subs({ksi: ksi_point, eta: eta_point}).evalf()
                        dN_i_deta = dN_deta_list[i].subs({ksi: ksi_point, eta: eta_point}).evalf()
                        dN_j_dksi = dN_dksi_list[j].subs({ksi: ksi_point, eta: eta_point}).evalf()
                        dN_j_deta = dN_deta_list[j].subs({ksi: ksi_point, eta: eta_point}).evalf()

                        integrand = (dN_i_dksi * dN_j_dksi + dN_i_deta * dN_j_deta) * J_inv
                        integral_value += integrand * gauss_weights[ksi_idx] * gauss_weights[eta_idx]

                K_local[i, j] = integral_value

        element_matrices.append(K_local)

    return element_matrices


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

def compute_jacobian(ksi, eta, x_coords, y_coords, dN_dksi_list, dN_deta_list):
    dx_dksi = sum(dN_dksi_list[j].subs({sp.symbols('ksi'): ksi, sp.symbols('eta'): eta}).evalf() * x_coords[j] for j in range(len(x_coords)))
    dx_deta = sum(dN_deta_list[j].subs({sp.symbols('ksi'): ksi, sp.symbols('eta'): eta}).evalf() * x_coords[j] for j in range(len(x_coords)))
    dy_dksi = sum(dN_dksi_list[j].subs({sp.symbols('ksi'): ksi, sp.symbols('eta'): eta}).evalf() * y_coords[j] for j in range(len(y_coords)))
    dy_deta = sum(dN_deta_list[j].subs({sp.symbols('ksi'): ksi, sp.symbols('eta'): eta}).evalf() * y_coords[j] for j in range(len(y_coords)))

    J_val = abs(dx_dksi * dy_deta - dx_deta * dy_dksi)
    return J_val

# Функція для обчислення частинних похідних
def compute_partial_derivatives(ap):
    ksi, eta = sp.symbols('ksi eta')

    num_nodes = 4 if ap == 1 else (9 if ap == 2 else 16)
    dN_dksi_list = []
    dN_deta_list = []

    for i in range(num_nodes):
        Ni = N(i, ksi, eta, ap)
        dN_dksi = sp.diff(Ni, ksi)
        dN_deta = sp.diff(Ni, eta)
        dN_dksi_list.append(dN_dksi)
        dN_deta_list.append(dN_deta)

    return dN_dksi_list, dN_deta_list

def apply_boundary_conditions(matrix, f_load, p, m, nodes, ug, ap = 1):
    bounds = get_boundary_points(p, m, ap)
    # print(bounds)
    n = len(ug)
    for i in range(n):
        for index in bounds[i]:
            for j in range(len(matrix[index])):
                matrix[index][j] = 0
            matrix[index][index] = 1
            f_load[index] = ug[i](nodes[index][0], nodes[index][1])
    return (matrix, f_load)

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

# Основна функція
def main():
    #vertices = [(0, 0), (1, 0.25), (1, 0.75), (0, 1)]
    vertices = [(0, 0), (1, 0), (1, 1), (0, 1)]
    element_type = 'D2QU4N'
    p = 5  # Кількість елементів по x
    m = 5  # Кількість елементів по y
    ap = 1
    ug = [ug_3, ug_2, ug_4, ug_1]
    NL, EL= uniform_mesh_with_vertices(vertices, p, m, element_type, ap)
    f_load = set_up_vector(f, NL, EL, p, m, ap)
    elem_matrices = compute_element_stiffness(EL, NL, ap)
    matrix =  assemble_global_stiffness_matrix(NL, EL, p, m, elem_matrices, ap)
    matrix, f_load = apply_boundary_conditions(matrix, f_load, p, m, NL, ug, ap)
    u = np.linalg.solve(matrix, f_load)

    m2d.plot_2d_solution(u, NL, EL)
    print("Координати вузлів:\n", NL)
    print("Елементи:\n", EL)
    print(f_load)
    print(matrix)
    print(u)


if __name__ == "__main__":
    main()
