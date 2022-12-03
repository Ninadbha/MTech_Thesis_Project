#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from sympy import *


# In[2]:


# coordinates of x
x = float(input("Enter the co-ordinate of x : "))
print("x :", x)


# In[3]:


# coordinates of position vector r
x_ = Symbol('x_')
# use x_ for making ux, uy and uz function instead of x
ux = float(input("Enter the directional displacement in  x direction : "))
uy = float(input("Enter the directional displacement in  y direction : "))
uz = float(input("Enter the directional displacement in  z direction : "))
r = np.array([[x+ux], [uy], [uz]])         # (1)
print("r :\n", r)


# In[4]:


# Jacobian
u_x = diff(ux, x_).subs({x_: x})
u_y = diff(uy, x_).subs({x_: x})
u_z = diff(uz, x_).subs({x_: x})

J = np.sqrt(float(np.square(1+u_x)) + float(np.square(u_y)) +float(np.square(u_z)))           # (2)
print("J :", J)


# In[5]:


# nodal displacements
ux1 = float(input("Enter the Nodal displacement of  x1 in local co-ordinate : "))
uy1 = float(input("Enter the Nodal displacement of  y1 in local co-ordinate : "))
uz1 = float(input("Enter the Nodal displacement of  z1 in local co-ordinate : "))
ux2 = float(input("Enter the Nodal displacement of  x2 in local co-ordinate : "))
uy2 = float(input("Enter the Nodal displacement of  y2 in local co-ordinate : "))
uz2 = float(input("Enter the Nodal displacement of  z2 in local co-ordinate : "))

ue = np.transpose(np.array([[ux1, uy1, uz1, ux2, uy2, uz2]]))         # (3)
print("ue :\n", ue)


# In[6]:


# section forces
fx1 = float(input("Enter the section forces in x1 direction : "))
fy1 = float(input("Enter the section forces in x1 direction : "))
fz1 = float(input("Enter the section forces in x1 direction : "))
fx2 = float(input("Enter the section forces in x1 direction : "))
fy2 = float(input("Enter the section forces in x1 direction : "))
fz2 = float(input("Enter the section forces in x1 direction : "))

fe = np.transpose(np.array([[fx1, fy1, fz1, fx2, fy2, fz2]]))           # (4)
print("fe :\n", fe)


# In[7]:


# Orthogonal Transformation matrix
T = np.array([[-1, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0],
              [0, 0, 0, 1, 0, 0],
              [0, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 1]])
print("T :\n", T)
print("Identity Matrix : \n", T @ np.transpose(T))


# In[8]:


# nodal displacements transformation to global coordinates
uge = np.transpose(T) @ ue         # (5)
print("uge :\n", uge)


# In[9]:


# sectional forces transformation to global coordinates
fge = np.transpose(T) @ fe               # (6)
print("fge :\n", fge)


# In[10]:


# uzg
# x_ = 2
l = float(input("enter the length : "))
uz1g = 5
uz2g = 6
uzg = (1 - x_/l)*uz1g + x*uz2g/l         # (8)
print("uzg :\n", uzg)


# In[11]:


# fwe
w = float(input("enter the gravity weight per unit length in N/m : "))
fwe = (J * w * l / 2 ) * np.transpose(np.array([[0, 0, 1, 0, 0, 1]]))            # (10)
print("fwe :\n", fwe)


# In[12]:


# pi_w
pi_w = np.linalg.det(np.transpose(-1 * uge) @ fwe)             # (9)
print("pi_w :", pi_w)


# In[13]:


# pi_e
pi_e = np.linalg.det(np.transpose(ue) @ fe)         # (11)
print("pi_e :", pi_e)


# In[14]:


s0 = float(input("enter the unstreched length(m) : "))
Tension = float(input("enter tension(N/m) : "))
EA = float(input("enter the value of EA : "))
pure_catenary = False        # other value can be True
s = s0                      # (13)
if not pure_catenary:
    s = s0 * (1 + Tension / EA)      # (14)

print("s :", s)


# In[15]:


A = (2 / (J+1)) * np.transpose([[-1, 0, 0, 1, 0, 0]])      # (17)
print("A :\n", A)


# In[16]:


# S = float(input("enter "))
# y = 14
# m = 15
C = (1 / ((J + 1) * l)) * np.array([[1, 0, 0, -1, 0, 0],        
                                    [0, 1, 0, 0, -1, 0],
                                    [0, 0, 1, 0, 0, -1],
                                    [-1, 0, 0, 1, 0, 0],
                                    [0, -1, 0, 0, 1, 0],
                                    [0, 0, -1, 0, 0, 1]])            # 18
print("C:\n", C)


# In[17]:


# pi_c
lmbda = 16
pi_c = lmbda * (np.linalg.det(np.transpose(ue) @ A) + np.linalg.det(np.transpose(ue) @ C @ ue) - (s-l))
print("pi_c:", pi_c)


# In[18]:


# pi
pi = (-1 * pi_w) + (-1 * pi_e) + pi_c            # (19)
print("pi :", pi)


# In[19]:


# del_pi_by_del_ue
Ke = lmbda * C            # (23)
B = A + C @ ue             # (24)
del_pi_by_del_ue = (Ke @ ue) + (T @ fwe) - (fe) + (B * lmbda)       # (22)
print("del_pi_by_del_ue :\n", del_pi_by_del_ue)


# In[20]:


# del_pi_by_del_lmbda
del_pi_by_del_lmbda = np.linalg.det((np.transpose(B) @ ue)) - (s-l)      # (25)
del_pi_by_del_lmbda


# In[21]:


Kge = np.transpose(T) @ Ke @ T
Bg = np.transpose(T) @ B
fge = np.transpose(T) @ fe
print("Kge :\n", Kge)
print("Bg :\n", Bg)
print("fge :\n", fge)


# In[22]:


K = np.r_[np.c_[Kge, Bg], np.c_[np.transpose(Bg), 0]]             # (35)
print("K :\n", K)
# len(K)  
# len(K[0])


# In[23]:


fn = np.r_[-1 * fwe, np.array([[s-l]])]              # (36)
print("fn :\n", fn)


# In[24]:


fge = Kge @ uge + Bg * lmbda + fwe                    # (37)
print("fge :\n", fge)


# In[ ]:





# In[ ]:




