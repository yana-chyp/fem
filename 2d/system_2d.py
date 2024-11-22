import nastia2d as n2d
import base_functions_2d as bf2d
import finite_element_2d as fe2d

def set_up_matrix(d1, d2, p, m, element_type='D2QU4N'):
    nodes, elements = n2d.uniform_mesh(d1, d2, p, m, element_type)
    #npe - nodes per element (assume all are the same)
    npe = len(elements[0])
    #let it be rectangle

    #assume nodes are equidistant
    h_x = d1/p
    h_y = d2/m
    J = h_x*h_y/4

    n = (p+1)*(m+1)
    matrix = [[0 for j in range(n)] for i in range(n)]

    em = fe2d.element_matrix()
    #ioe - index of element
    for ioe in range(len(elements)):
        #noe - nodes of element
        noe = elements[ioe]
        # print(noe)

        # im = noe[0]//(p+1)
        # jm = noe[0] - im * (p+1)
        for i in range(len(noe)):
            for j in range(len(noe)):
                im = noe[i]
                jm = noe[j]
                value = em[i][j]*J

                # print("(", i, ",  ", j, ") -> (", im, ", ", jm, "); value = ", f"{value:.4f}")
                matrix[im][jm] += value
            # print('[' + ', '.join([f"{el:.8f}" for el in matrix[i]]) + ']')
    return matrix




def get_boundary_points(p, m):
    bounds = []
    for j in range(p+1):
        bounds.append(j)
    for i in range(1, m+1):
        index = (i+1)*(p+1)-1
        bounds.append(index)
    for j in reversed(range(p)):
        index = p*(m+1) - 1 + j
        bounds.append(index)
    for i in reversed(range(1, m)):
        index = i*(p+1)
        bounds.append(index)
    return bounds




#let ug_i = f(x) = y - function defining the boundary conditions
#here i = 1,...,4 means 4 intervals on the plane intersecting the Omega area
# def get_boundary_points(NL , ug_1, ug_2, ug_3, ug_4):
    # gamma_1 = []
    # gamma_2 = []
    # gamma_3 = []
    # gamma_4 = []
    # for i in range(len(NL)):
    #     if NL[i][1]==ug_1(NL[i][0]):
    #         gamma_1.append(NL[i])
    #     elif NL[i][1]==ug_2(NL[i][0]):
    #         gamma_2.append(NL[i])
    #     elif NL[i][1]==ug_3(NL[i][0]):
    #         gamma_3.append(NL[i])
    #     elif NL[i][1]==ug_4(NL[i][0]):
    #         gamma_4.append(NL[i])

    # return[gamma_1, gamma_2, gamma_3, gamma_4]



