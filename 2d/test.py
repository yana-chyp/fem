import numpy as np
import sympy as sp
import mesh_2d as m2d
import system_2d as s2d
import finite_element_2d as f2e
import graph2d as g2d
from enum import Enum
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def ug_1(x, y):
    return 0
    # return y**2
    # return y
    # return y**2


def ug_2(x, y):
    return 0
    # return y**2
    # return y+1
    # return y**2 + 1


def ug_3(x, y):
    return 0
    # return x**2
    # return x


def ug_4(x, y):
    return 0
    # return x**2+1
    # return x+1


def k1(x, y):
    return 1


def k2(x, y):
    return 1


# test1
def ft1(x, y, A=100, x0=0.5, y0=0.5, const=1):
    return A * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / const)


def ft2(x, y, A=100, x0=0.5, y0=0.5, const=0.1):
    return A * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / const)


# test2
sourcest3 = [(0.5, 0.5)]
strengthst3 = [100]

sourcest4 = [(0.5, 0.5), (0.2, 0.8)]
strengthst4 = [100, 50]

atol = 0.05

ft3 = lambda x, y: sum(
    s for (x0, y0), s in zip(sourcest3, strengthst3)
    if np.isclose(x, x0, atol=atol) and np.isclose(y, y0, atol=atol)
)

ft4 = lambda x, y: sum(
    s for (x0, y0), s in zip(sourcest4, strengthst4)
    if np.isclose(x, x0, atol=atol) and np.isclose(y, y0, atol=atol)
)


# test3
def k1t5(x, y):
    return 0.1 if 0.45 < x < 1.55 else 10.0


def k2t5(x, y):
    return 1.0


def ft5(x, y):
    return 100 * np.exp(-((x - 0.5) ** 2 + (y - 0.5) ** 2) / 0.001)


def k1t6(x, y):
    if 0.45 < x < 1.55:
        return 0.05
    return 1.0


def k2t6(x, y):
    if 0.65 < y < 1.75:
        return 0.01
    return 1.0


def ft6(x, y):
    if x < 0.3 and y < 0.3:
        return 20
    elif x > 0.7 and y > 0.7:
        return 5
    else:
        return 0


# ver
# def fv(x, y):
#     return 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
def fv(x, y):
    return 2 * sp.pi ** 2 * sp.sin(sp.pi * x) * sp.sin(sp.pi * y)


def exact_solution(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y)


def main():
    verticest1 = [(0, 0), (1.5, 0), (1, 1), (0, 0.75)]
    verticest2 = [(0, 0), (1.5, 0), (1, 1), (0, 0.75)]
    verticest3 = [(0, 0), (1, 0), (1.1, 1), (0, 0.9)]
    verticest4 = [(0, 0), (1, 0), (1.1, 1), (0, 0.9)]
    verticest5 = [(0, 0), (1.5, 0), (1, 1), (0, 0.75)]
    verticest6 = [(0, 0), (1, 0), (1, 1), (0, 1)]
    verticesv = [(0, 0), (1, 0), (1, 1), (0, 1)]

    ap = 1

    ug = [
        [f2e.TypeOfBoundCond.DIRICHLET, ug_3],  # нижнє (y = 0)
        [f2e.TypeOfBoundCond.DIRICHLET, ug_2],  # ліве (x = 0)
        [f2e.TypeOfBoundCond.DIRICHLET, ug_4],  # праве (x = 1)
        [f2e.TypeOfBoundCond.DIRICHLET, ug_1]  # верхнє (y = 1)
    ]

    element_type = 'D2QU4N'
    p = 5
    m = 5


if __name__ == '__main__':
    main()
