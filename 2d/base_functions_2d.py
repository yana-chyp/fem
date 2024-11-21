import sympy as sp


ksi_left = -1
ksi_right = 1
eta_left = -1
eta_right = 1

def get_base_functions(m = 1):
    #for now gives correct results only for m = 1
    #greater values of m need work
    n = (m+1)*(m+1)     #number of nodes
    base_functions = [0 for k in range(n)]
    ksi = sp.symbols('ksi')
    eta = sp.symbols('eta')

    h_ksi = (ksi_right - ksi_left)/m
    h_eta = (eta_right - eta_left)/m

    submesh = [[ksi_left + k*h_ksi, eta_left + l*h_eta] for l in range(m+1) for k in range(m+1)]
    # print(submesh)
    for i in range(n):
        expression = 1
        for k in range(m+1):
            if submesh[k][0] != submesh[i][0]:
                expression *= ( (ksi - submesh[k][0])/(submesh[i][0]-submesh[k][0]) )
        for l in range(m+1):
            if  submesh[l*(m+1)][1] != submesh[i][1]:
                expression *= ( (eta - submesh[l*(m+1)][1])/(submesh[i][1]-submesh[l*(m+1)][1]))


        base_functions[i] = sp.simplify(expression)
    #     for j in range(m+1):
    #         expression = 1
    #         for k in range(m + 1):
    #             if k != i:
    #                 expression *= (ksi - submesh[k][j][0]) / (submesh[i][j][0] - submesh[k][j][0])
    #         for l in range(m + 1):
    #             if l != j:
    #                 expression *= (eta - submesh[i][l][1]) / (submesh[i][j][1] - submesh[i][l][1])
    #         base_functions[i][j] = sp.simplify(expression)



    return base_functions

