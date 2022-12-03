#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


def func(H):
    return (((w*xL)/(np.arccosh(1 - (w*w)*(zL*zL-st*st)/(2*H*H)))) - H)

def derivFunc(H):
    return (((-1*w*xL)/(np.square(np.arccosh(1 - (w*w)*(zL*zL-st*st)/(2*H*H))))) * (1/(np.sqrt(np.square(1 - (w*w)*(zL*zL-st*st)/(2*H*H))-1))) * ((w*w)*(zL*zL-st*st)/(H*H*H)) - 1)

def newtonRaphson(H):
    h = func(H) / derivFunc(H)
    while abs(h) >= 0.001:
        h = func(H)/derivFunc(H)
        H = H - h
    return H


# In[4]:


theta = int(input("enter theta in degrees : "))  # in degree
l = int(input("enter the length (m) : "))
xL = l * np.cos(theta * np.pi / 180)  # when theta is in degree
zL = l * np.sin(theta * np.pi / 180)  # when theta is in degree
print(zL)
w = int(input("enter the weight per unit length(N/m) : "))
st = int(input("enter the total unstreched length(m)  : "))
H = newtonRaphson(-20)    # initial guess -20
xs = (H/w) * np.log((w*(zL-st))/(H * (np.exp(-1*w*xL/H)-1)))
print(f"xs : {xs} m")


# In[4]:


x = np.linspace(0, 1000, 50)
z = np.array([])
for a in x:
    b = (H/w)*(np.cosh((w/H)*(xs-a))) - (H/w)*(np.cosh(w*xs/H))
    z = np.append(z, b)


# In[5]:


plt.figure(figsize=(5, 5))
plt.scatter(x, z, label=" Analytic  ")
plt.plot(x, z,'ro')
plt.xlabel("x (m)")
plt.ylabel("z (m)")
plt.title(f"Analytic results for equilibrium \n position of example catenary(theta = {theta} degrees)")
plt.legend()
plt.show()


# In[6]:


x = np.linspace(0, 1000, 50)
t = np.array([])
for a in x:
    v = -1 * H * np.sinh((w/H)*(xs-a))
    b = np.sqrt(np.square(H) + np.square(v))
    t = np.append(t, b)


# In[7]:


plt.figure(figsize=(5, 5))
plt.scatter(x, t, label=" Analytic  ")
plt.plot(x, t ,'go')
plt.xlabel("x (m)")
plt.ylabel("T (N)")
plt.title(f"Analytic results for tension distribution of \n  example catenary(theta = {theta} degrees)")
plt.legend()
plt.show()


# In[ ]:




