import sympy as sp


ksi_left = -1
ksi_right = 1

def get_x(ksi, i, mesh):
    return (ksi+1)/2 * (mesh[i+1] - mesh[i]) + mesh[i]

def get_ksi(x, i, mesh):
    return 2*(x - mesh[i])/(mesh[i+1] - mesh[i]) - 1

def get_base_functions(m):
    #here i is the index of the left end of the interval
    base_functions = []
    ksi = sp.symbols('ksi')



    h = (ksi_right - ksi_left)/m
    submesh = [ksi_left + k*h for k in range(m+1)]
    # print(i, submesh)
    for k in range(m+1):
        expression = 1
        for j in range(m+1):
            if k!=j:
                expression *= (ksi - submesh[j])/(submesh[k] - submesh[j])
        base_functions.append(expression)



    #func = sp.lambdify(x, expression)  ?
    return base_functions
