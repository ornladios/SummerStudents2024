import sympy as sp

import numpy as np

 

# Define the variables

x, y, z = sp.symbols('x y z')

 
x_prime = (x+ts )
y_prime =
z_prime = 
# Define the function

f = sp.exp((x**2 + y**2 + z**2)/(8*(np.pi**3)))

 

# Compute the Laplacian

laplacian = sp.diff(f, x, x) + sp.diff(f, y, y) + sp.diff(f, z, z)

 

# Print the result

print(laplacian)