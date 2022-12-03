#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import cmath
import sympy as smp
from sympy import*


# In[2]:


def func(H):
    return (((w*xL)/(np.arccosh(1 - (w*w)*(zL*zL-st*st)/(2*H*H)))) - H)

def derivFunc(H):
    return (((-1*w*xL)/(np.square(np.arccosh(1 - (w*w)*(zL*zL-st*st)/(2*H*H))))) * (1/(np.sqrt(np.square(1 - (w*w)*(zL*zL-st*st)/(2*H*H))-1))) * ((w*w)*(zL*zL-st*st)/(H*H*H)) - 1)


# In[3]:


w=int(input("enter the weight in N\m : "))

L = int(input("enter the length : "))
print(f"length is {L}")

ϴ = int(input("enter the angle in degrees : "))
# ϴ = ϴ *(180/math.pi)
# ϴ = ϴ *(math.pi/180)
print(f"theta is {ϴ}")
# xL= int(input("enter the length between the"))

xL = L*math.cos(ϴ)
print(f"xL is {xL}")

zL = L*math.sin(ϴ)
print(zL)
# zL=3

st = int(input("Enter the Total unstretched length: "))

def newtonRaphson(H):
    h = func(H) / derivFunc(H)
    while abs(h) >= 0.001:
        h = func(H)/derivFunc(H)
        H = H - h
        print(abs(h))
    print("The value of the root is : ",
                             "%.4f"% H)

H0 = -10 # Initial values assumed
newtonRaphson(H0)


# In[4]:


# func(1283520.1392)
H = func( -719383.7793)
print(f"H is : {H}")


# In[ ]:





# In[ ]:




