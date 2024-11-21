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
    J_x = 2/h_x
    J_y = 2/h_y

    matrix = [[0 for j in range(m+1)] for i in range(p+1)]

    em = fe2d.element_matrix()
    #ioe - index of element
    for ioe in range(len(elements)):
        #noe - nodes of element
        noe = elements[ioe]
        # print(noe)
        for i in range(len(noe)):
            for j in range(len(noe)):
                im = noe[0]//(p+1)
                jm = noe[0] - im * (p+1)
                value = em[i][j][0]*J_x**2 + em[i][j][1]*J_y**2
                print("i = ", i, " j = ", j, " value = ", value)
                matrix[im][jm] += value

    return matrix