import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def fourth_derivative(f, axis, h):
    h =2*np.pi/29
    #print(h)
    d2f = np.gradient(np.gradient(f, h, axis=axis), h, axis=axis)
    d4f = np.gradient(np.gradient(d2f, h, axis=axis), h, axis=axis)
    return d4f

h = 2*np.pi/29
#print(h)
smin, smax = 0, 2*np.pi
npt = int((smax-smin)/h)+1
print(npt, npt**3 * 8 / 1024/1024/1024)
x = np.linspace(smin, smax, npt)
y = np.linspace(smin, smax, npt)
z = np.linspace(smin, smax, npt)
xs, ys,zs = np.meshgrid(x,y,z)
scalar = (2*np.pi)**3
f = np.exp((xs**2+ys**2+zs**2)/scalar)
L_f = (6.50100920809908e-5*xs**2 + 0.00806288360829987)*np.exp(0.00403144180414994*xs**2 + 0.00403144180414994*ys**2 + 0.00403144180414994*zs**2) + (6.50100920809908e-5*ys**2 + 0.00806288360829987)*np.exp(0.00403144180414994*xs**2 + 0.00403144180414994*ys**2 + 0.00403144180414994*zs**2) + (6.50100920809908e-5*zs**2 + 0.00806288360829987)*np.exp(0.00403144180414994*xs**2 + 0.00403144180414994*ys**2 + 0.00403144180414994*zs**2)
# truncation error = L_F - numerical solction from CalcLaplsian 
#fig, axs = plt.subplots(1)
#axs[0].imshow(f[20,:,:])
#plt.imshow(L_f[:,:,20])
#plt.colorbar()
#plt.savefig('test_Lap.png')

import adios2 
import sys
from adios2 import Stream 

fname = sys.argv[1]
var = sys.argv[2]

#T_err = h**2/12 *( fiv(x) + fiv(y) + fiv(z))
num = (h**2)/12
fx4 = fourth_derivative(f, axis=0, h=h)
fy4 = fourth_derivative(f, axis=1, h=h)
fz4 = fourth_derivative(f, axis=2, h=h)
T_err = num*(fx4 + fy4 + fz4)

with Stream(fname, 'r') as file:
    for _ in file.steps():
        step = file.current_step()
        data = file.read(var)
        errors = np.abs(L_f - data)
       
        errors[:,  :, 0] = 0
        errors[:,  :, -1] = 0
        errors[:,  0, :] = 0
        errors[:, -1, :] = 0
        errors[0,  :, :] = 0
        errors[-1,  :, :] = 0
   
        print("T_err from data calulated: ", errors.max())
        print("T_err from python script ", np.max(T_err))


# read Laplace in
# remove 0's from Laplace
# subtract numeric vals - Laplace read vals
# L_F - Laplace read in