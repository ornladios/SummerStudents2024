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

        rmse = np.sqrt(np.mean(errors**2))
        print(f"rmse: {rmse}")
        nomralized = rmse/(data.max()-data.min())
        print("T_err from data calulated: ", errors.max())
        print("T_err from python script ", np.max(T_err))
        print(f"normalized error: {nomralized}")


# read Laplace in
# remove 0's from Laplace
# subtract numeric vals - Laplace read vals
# L_F - Laplace read in

# error bound 1*10^-4 
# qunaitzation error normailzed, second is compressed error, use same code to calculate the same error 
# read compressed data 
# write is uncormpressed 
# when you read in python if the code is correct when adios reads the compress it combines the read and read reconstered, 

# 1. C++ --> write out mgard.bp -data f, L_f, --> size 30^3 use FIRST function 
# 2. xml file define operatpor (mgard), error bound (0.0001) 0r (0.00006)
# compress L_f 
# measerure the error in L_f small size
# 3. python  script: read adios -mgard, adios will return rct_data same size as orginal (check) 
# 4. compute comp_err = np.sqrt(np.mean(diff**2)) , where diff = rct_L - L_f 
#5 check if comp comp_err < T_err, or comp_err/norm < T_err/norm. where norm =L.max()-min()
# 6. repeat across multiple time steps 
# Fun = exp(x^2+y^2+z^2)/scalar
#change x,y,z to x',y',and z' : x' =(x+ts*delta),delta =1.5 --> 20 timespets. 
#
# 7. plot: x-axis = ts, plot both T_err(ts) and comp_err(ts)

"""
1. C++ to write out mgard.bp – f, L_f   30x30x30
                 - xml file to define operator (mgard), error bound  change to 0.0001 or 0.00006
                 - compress L_f, f ??? small size
2. Python script: file.read adios – mgard.bp, adios will return rct_L, rct_f  30x30x30
3. comp_err = np.sqrt(np.mean(diff**2)), where diff = rct_L – L_f
4. check if comp_err < T_err, or if comp_err/norm < T_err / norm, where norm = L.max()-L.min()
5. Repeat this for multiple timesteps
Fun = exp((x**2 + y**2 + z**2) / scalar)
Change x, y, z to x’, y’, and z’: x’ = (x+ts*delta), delta = 1.5  20 timesteps , same thing for y’ and z’
Python script: change fun to calculate the analytical solution
6. Plot: x-axis = ts, plot both T_err(ts) and comp_err(ts)
7. compression ratio: averave compression ration over ts
write out L_f --> mgard.bp 
compression ratio = 30^3*8/size



Backup taks:
vis 1-2 slice two images
add a plot: compresssion ratio vs error

motivate 2 things
ensure that this is an error bounded compresion 
we can signaiffecntly compression this and still have accurate data 

add compression ratio!!!!!!!!!!!!!!!!
what is benift of using 

eb = 0.01
x-axis = [0.01, 0.005, 0.001, 0.0005, 0.00001]
y-axis = compressisoin ratio 

plotimg compresison aratio vs erro bound 


"""
