import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import dblquad
import scipy.integrate as scin
from mpl_toolkits.mplot3d import Axes3D


ksi_left = -1
ksi_right = 1
eta_left = -1
eta_right = 1

def f(x, y):
    return (1 + x**2)*(1 + 2*y**2)

def uniform_mesh(d1, d2, p, m, element_type):
    PD = 2  # Простір (x, y)
    q = np.array([[0, 0], [d1, 0], [0, d2], [d1, d2]])  # 4 кути прямокутника

    NoN = (p + 1) * (m + 1)  # Кількість вузлів
    NoE = p * m  # Кількість елементів (для прямокутників)
    NPE = 4 if element_type == 'D2QU4N' else 3  # Кількість вузлів на елемент

    NL = np.zeros([NoN, PD])
    a = (q[1, 0] - q[0, 0]) / p  # Інкременти по x
    b = (q[2, 1] - q[0, 1]) / m  # Інкременти по y

    n = 0
    for i in range(m + 1):
        for j in range(p + 1):
            NL[n, 0] = q[0, 0] + j * a  # Значення x
            NL[n, 1] = q[0, 1] + i * b  # Значення y
            n += 1

    if element_type == 'D2QU4N':  # Прямокутники
        EL = np.zeros([NoE, NPE], dtype=int)
        for i in range(m):
            for j in range(p):
                n1 = i * (p + 1) + j
                n2 = n1 + 1
                n3 = n1 + (p + 1)
                n4 = n3 + 1
                EL[i * p + j] = [n1, n2, n4, n3]
    elif element_type == 'D2TR3N':  # Трикутники
        EL = np.zeros([NoE * 2, NPE], dtype=int)  # Подвійна кількість елементів для трикутників
        e = 0
        for i in range(m):
            for j in range(p):
                n1 = i * (p + 1) + j
                n2 = n1 + 1
                n3 = n1 + (p + 1)
                n4 = n3 + 1
                # Перший трикутник
                EL[e] = [n1, n2, n4]
                e += 1
                # Другий трикутник
                EL[e] = [n1, n4, n3]
                e += 1
    else:
        raise ValueError("Неправильний тип елемента. Використовуйте 'D2QU4N' або 'D2TR3N'.")

    return NL, EL

def get_base_functions(m=1):
    # Повертає базисні функції для m+1 x m+1 вузлів
    n = (m + 1) * (m + 1)  # Кількість вузлів
    base_functions = [0 for _ in range(n)]
    ksi, eta = sp.symbols('ksi eta')

    # Координати вузлів у локальній системі
    h_ksi = (ksi_right - ksi_left) / m
    h_eta = (eta_right - eta_left) / m

    submesh = [[ksi_left + j * h_ksi, eta_left + i * h_eta]
               for i in range(m + 1) for j in range(m + 1)]

    for i in range(n):
        expression = 1
        for k in range(n):
            if k != i:
                dx = submesh[i][0] - submesh[k][0]
                dy = submesh[i][1] - submesh[k][1]
                if dx != 0:
                    expression *= (ksi - submesh[k][0]) / dx
                if dy != 0:
                    expression *= (eta - submesh[k][1]) / dy
        base_functions[i] = expression  # Без спрощення на цьому етапі

    return base_functions


def N(i, ksi, eta):
    """Базисні функції для 4-кутних елементів."""
    if i == 0:
        return (1 - ksi) * (1 - eta) / 4
    elif i == 1:
        return (1 + ksi) * (1 - eta) / 4
    elif i == 2:
        return (1 + ksi) * (1 + eta) / 4
    elif i == 3:
        return (1 - ksi) * (1 + eta) / 4
    return 0

def calculate_b(NL, EL, p, m):
    """Обчислення вектора правої частини для заданої сітки."""
    b = np.zeros(NL.shape[0])  # Ініціалізація глобального вектора правої частини

    # Інтегрування для кожного елемента
    for e in range(EL.shape[0]):
        nodes = EL[e]  # Вузли елемента
        x_coords = NL[nodes, 0]
        y_coords = NL[nodes, 1]

        # Вирази для перетворення координат
        def x_expr(ksi, eta):
            return sum(N(i, ksi, eta) * x_coords[i] for i in range(4))

        def y_expr(ksi, eta):
            return sum(N(i, ksi, eta) * y_coords[i] for i in range(4))

        # Детермінант Якобіана
        def det_jacobian(ksi, eta):
            J11 = sum((1 - eta if i % 2 == 0 else 1 + eta) * x_coords[i] / 4 for i in range(4))
            J12 = sum((1 - ksi if i < 2 else 1 + ksi) * x_coords[i] / 4 for i in range(4))
            J21 = sum((1 - eta if i % 2 == 0 else 1 + eta) * y_coords[i] / 4 for i in range(4))
            J22 = sum((1 - ksi if i < 2 else 1 + ksi) * y_coords[i] / 4 for i in range(4))
            return abs(J11 * J22 - J12 * J21)
        # def det_jacobian(ksi, eta, x_coords, y_coords):
        #     """Обчислення детермінанта Якобіана."""
        #     # Похідні базисних функцій
        #     dN_dksi = [sp.diff(N(i, sp.Symbol('ksi'), sp.Symbol('eta')), 'ksi') for i in range(4)]
        #     dN_deta = [sp.diff(N(i, sp.Symbol('ksi'), sp.Symbol('eta')), 'eta') for i in range(4)]
        #
        #     # Елементи Якобіана
        #     J11 = sum(dN_dksi[i].subs({'ksi': ksi, 'eta': eta}) * x_coords[i] for i in range(4))
        #     J12 = sum(dN_deta[i].subs({'ksi': ksi, 'eta': eta}) * x_coords[i] for i in range(4))
        #     J21 = sum(dN_dksi[i].subs({'ksi': ksi, 'eta': eta}) * y_coords[i] for i in range(4))
        #     J22 = sum(dN_deta[i].subs({'ksi': ksi, 'eta': eta}) * y_coords[i] for i in range(4))
        #
        #     # Детермінант Якобіана
        #     return abs(J11 * J22 - J12 * J21)

        # Інтегральна функція для вузла i
        def integrand(ksi, eta, i):
            x = x_expr(ksi, eta)
            y = y_expr(ksi, eta)
            return N(i, ksi, eta) * f(x, y) * det_jacobian(ksi, eta)
        # def integrand(ksi, eta, i, x_coords, y_coords):
        #     """Інтегральна функція для вузла i."""
        #     x = sum(N(j, ksi, eta) * x_coords[j] for j in range(4))
        #     y = sum(N(j, ksi, eta) * y_coords[j] for j in range(4))
        #     det_J = det_jacobian(ksi, eta, x_coords, y_coords)
        #     return N(i, ksi, eta) * f(x, y) * det_J

        # Обчислення локальних інтегралів для кожного вузла елемента
        for i in range(4):
            integral, _ = dblquad(
                lambda eta, ksi: integrand(ksi, eta, i),
                -1, 1,  # Межі для ksi
                lambda ksi: -1, lambda ksi: 1  # Межі для eta
            )
            b[nodes[i]] += integral

    return b




def plot_2d_solution(u, NL, EL, exact_solution=None):
    """
    Побудова графіка наближеного розв'язку в 2D-просторі.

    u - знайдений вектор розв'язку (значення на вузлах).
    NL - координати вузлів (матриця розмірності [NoN x 2]).
    EL - елементи сітки (матриця розмірності [NoE x NPE]).
    exact_solution - функція для точного розв'язку (опціонально).
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Нанесення вузлів
    x = NL[:, 0]  # x-координати вузлів
    y = NL[:, 1]  # y-координати вузлів
    z = u  # Значення наближеного розв'язку в вузлах

    # Побудова триангуляції (за елементами EL)
    ax.plot_trisurf(x, y, z, cmap='viridis', alpha=0.8, edgecolor='gray', label='Наближений розв\'язок')

    # Якщо заданий точний розв'язок
    if exact_solution:
        x_exact = np.linspace(np.min(x), np.max(x), 50)
        y_exact = np.linspace(np.min(y), np.max(y), 50)
        X, Y = np.meshgrid(x_exact, y_exact)
        Z = exact_solution(X, Y)
        ax.plot_surface(X, Y, Z, cmap='plasma', alpha=0.4, label='Точний розв\'язок')

    # Налаштування графіка
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u(x, y)')
    ax.set_title('Наближений розв\'язок у 2D-просторі')
    plt.legend()
    plt.show()

def uniform_mesh_level1(d1, d2, p, m, element_type):
    PD = 2  # Простір (x, y)
    q = np.array([[0, 0], [d1, 0], [0, d2], [d1, d2]])  # 4 кути прямокутника

    NoN = (p + 1) * (m + 1)  # Кількість вузлів
    NoE = p * m  # Кількість елементів (для прямокутників)
    NPE = 4 if element_type == 'D2QU4N' else 3  # Кількість вузлів на елемент

    NL = np.zeros([NoN, PD])
    a = (q[1, 0] - q[0, 0]) / p  # Інкременти по x
    b = (q[2, 1] - q[0, 1]) / m  # Інкременти по y

    n = 0
    for i in range(m + 1):
        for j in range(p + 1):
            NL[n, 0] = q[0, 0] + j * a  # Значення x
            NL[n, 1] = q[0, 1] + i * b  # Значення y
            n += 1

    EL = np.zeros([NoE, NPE], dtype=int)
    for i in range(m):
        for j in range(p):
            n1 = i * (p + 1) + j
            n2 = n1 + 1
            n3 = n1 + (p + 1)
            n4 = n3 + 1
            EL[i * p + j] = [n1, n2, n3, n4]

    return NL, EL

def uniform_mesh_level2(d1, d2, p, m, element_type):
    PD = 2  # Простір (x, y)
    nodes_x = 2 * p + 1  # Подвоєна кількість вузлів по x (основні + середні)
    nodes_y = 2 * m + 1  # Подвоєна кількість вузлів по y (основні + середні)

    NL = np.zeros([nodes_x * nodes_y, PD])
    x_coords = np.linspace(0, d1, nodes_x)  # Координати вузлів по x
    y_coords = np.linspace(0, d2, nodes_y)  # Координати вузлів по y

    n = 0
    for y in y_coords:
        for x in x_coords:
            NL[n, :] = [x, y]
            n += 1

    NoE = p * m  # Кількість елементів
    EL = np.zeros([NoE, 9], dtype=int)  # 8 вузлів на елемент (лише на краях)

    e = 0
    for i in range(m):
        for j in range(p):
            n1 = i * 2 * nodes_x + j * 2
            n2 = n1 + 2
            n3 = n1 + 2 * nodes_x
            n4 = n3 + 2
            n5 = n1 + 1
            n6 = n3 + 1
            n7 = n1 + nodes_x
            n8 = n7 + 2
            n9 = n7 + 1
            # [n1, n5, n2, n8, n4, n6, n3, n7, n9]
            EL[e, :] = [n1, n5, n2, n7, n9, n8, n3, n6, n4]
            e += 1

    return NL, EL

def uniform_mesh_level3(d1, d2, p, m, element_type):
    PD = 2  # Простір (x, y)
    nodes_x = 3 * p + 1  # Потроєна кількість вузлів по x
    nodes_y = 3 * m + 1  # Потроєна кількість вузлів по y

    NL = np.zeros([nodes_x * nodes_y, PD])
    x_coords = np.linspace(0, d1, nodes_x)  # Координати вузлів по x
    y_coords = np.linspace(0, d2, nodes_y)  # Координати вузлів по y

    n = 0
    for y in y_coords:
        for x in x_coords:
            NL[n, :] = [x, y]
            n += 1

    NoE = p * m  # Кількість елементів
    EL = np.zeros([NoE, 12], dtype=int)  # 12 вузлів на елемент

    e = 0
    for i in range(m):
        for j in range(p):
            n1 = i * 3 * nodes_x + j * 3
            n2 = n1 + 3
            n3 = n1 + 3 * nodes_x
            n4 = n3 + 3
            n5 = n1 + 1
            n6 = n1 + 2
            n7 = n3 + 1
            n8 = n3 + 2
            n9 = n1 + nodes_x
            n10 = n9 + 1
            n11 = n9 + 2
            n12 = n9 + 3
            EL[e, :] = [n1, n5, n6, n2, n12, n8, n7, n11, n3, n9, n10, n4]
            e += 1

    return NL, EL
def visualize_mesh(NL, EL, element_type):
    NoN = NL.shape[0]  # Кількість вузлів
    NoE = EL.shape[0]  # Кількість елементів

    plt.figure(figsize=(10, 10))

    # Анотація вузлів
    for i in range(NoN):
        plt.scatter(NL[i, 0], NL[i, 1], color='black', s=10)  # Вузли
        plt.text(NL[i, 0], NL[i, 1], str(i + 1), color='red', fontsize=8)

    # Побудова елементів
    for j in range(NoE):
        nodes = EL[j]
        x = NL[nodes, 0]
        y = NL[nodes, 1]
        x = np.append(x, x[0])  # Замкнути контур
        y = np.append(y, y[0])
        plt.plot(x, y, color='blue')
        # Анотація елемента
        cx = np.mean(x[:-1])
        cy = np.mean(y[:-1])
        plt.text(cx, cy, str(j + 1), color='green', fontsize=10)

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'Сітка: {element_type}')
    plt.grid()
    plt.axis('equal')
    plt.show()


def main():
    d1 = 1  # Довжина по x
    d2 = 1  # Довжина по y
    p = 2  # Поділ по x
    m = 2  # Поділ по y

    # Рівень 1
    NL1, EL1 = uniform_mesh_level1(d1, d2, p, m, 'D2QU4N')
    visualize_mesh(NL1, EL1, 'D2QU4N')

    # Рівень 2
    NL2, EL2 = uniform_mesh_level2(d1, d2, p, m, 'D2QU4N')
    # visualize_mesh(NL2, EL2, 'D2QU4N')

    # Рівень 3
    NL3, EL3 = uniform_mesh_level3(d1, d2, p, m, 'D2QU4N')
    # visualize_mesh(NL3, EL3, 'D2QU4N')


if __name__ == "__main__":
    main()