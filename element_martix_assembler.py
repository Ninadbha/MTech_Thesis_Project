#!/usr/bin/env python
# coding: utf-8

# In[1]:


# 
import numpy as np

'''type_test = 1
elements_test = 2
A_test = "200,200"
E_test = "200000.0,200000"
L_test = "200.0,200"
Q_test = "0,x,x"
F_test = "0,0,100.0"

type_test = 0
elements_test = 2
A_test = 200
E_test = 200000.0
L_test = 200
Q_test = 0
F_test = "0,0,90000"'''

def solve(type, elements, A, E, L, Q, F):

    if type == 0:

        L = L/elements
        K = A*E/L

        Global_Q = np.zeros((1, elements+1)).flatten()
        Global_Q = Global_Q.astype(str)

        if float(Q) == 0:
            for i in range(elements):
                Global_Q[i+1] = 'x'
        else:
            for i in range(elements-1):
                Global_Q[i+1] = 'x'
                Global_Q[elements] = str(Q)

        #Global_Q = Global_Q.reshape(elements + 1, 1)

        Global_F = np.array(F.split(","))
        Global_F = Global_F.reshape(elements + 1, 1).astype(float)


        # Assembling the global stiffness matrix

        Global_K = np.zeros((elements + 1, elements + 1))

        for i in range(elements):
            for row in range(2):
                for column in range(2):
                    if row == column:
                        Global_K[row + i][column + i] += K
                    else:
                        Global_K[row + i][column + i] -= K

    else:
        A = np.array(A.split(","))
        A = A.astype(float)

        E = np.array(E.split(","))
        E = E.astype(float)

        L = np.array(L.split(","))
        L = L.astype(float)

        Global_Q = np.array(Q.split(","))
        print(Global_Q)
        #Global_Q = Q.reshape(elements + 1, 1)

        Global_F = np.array(F.split(","))
        Global_F = Global_F.reshape(elements + 1, 1).astype(float)


        K = np.array([])

        for a, e, l in zip(A, E, L):
            K = np.append(K, a*e/l)

        #Assembling the global stiffness matrix

        Global_K = np.zeros((elements+1, elements+1))

        for i in range(elements):
            for row in range(2):
                for column in range(2):
                    if row == column:
                        Global_K[row+i][column+i] += K[i]
                    else:
                        Global_K[row + i][column + i] -= K[i]

#--------------------------matrix calculations-------------------------#

    print("\nThe global stiffness matrix is: ")
    print((Global_K))

    X = np.where(Global_Q != 'x')
    rows = np.unique(X[0]).astype(int)


    # Global force matrix

    print("\nThe global force matrix is:")
    print((Global_F))

    for i in rows:
        for f in range(elements+1):
            Global_F[f] -= Global_K[f][i]*float(Global_Q.reshape(elements+1, 1)[i])

    F = Global_F

    F = np.delete(F, rows, 0)

    print("\nThe global force matrix after elimination is:")
    print((F))

    Global_K = np.delete(Global_K, rows, 0)
    Global_K = np.delete(Global_K, rows, 1)

    print("\nThe global stiffness matrix after elimination is: ")
    print((Global_K))

    unknown_disp = np.linalg.solve(Global_K,F)

    i = 0

    for row in range(elements+1):
        if row in rows:
            continue
        else:
            Global_Q[row] = str(0)
            i += 1

    Global_Q = Global_Q.astype(float)

    i = 0

    for row in range(elements+1):
        if row in rows:
            continue
        else:
            Global_Q[row] = unknown_disp[i]
            i += 1

    print("\n The global nodal displacement matrix with all values is:")
    print(Global_Q)

    stress = np.array([])
    strain = np.array([])

    if type == 0:
        for element in range(elements):
            stress = np.append(stress, E * ((Global_Q[element + 1] - Global_Q[element]) / L))
            strain = np.append(strain, (Global_Q[element + 1] - Global_Q[element]) / L)

    else:
        for element in range(elements):
            stress = np.append(stress, E[element]*((Global_Q[element+1]-Global_Q[element])/L[element]))
            strain = np.append(strain, (Global_Q[element + 1] - Global_Q[element]) / L[element])

    print("\nThe stress matrix is:")
    print(stress)
    print("\nThe strain matrix is:")
    print(strain)

    return Global_Q.flatten(), Global_F.reshape(1, elements+1).flatten(), stress, strain

#solve(type_test, elements_test, A_test, E_test, L_test, Q_test, F_test)

