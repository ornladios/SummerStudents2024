import sys
import numpy as np
from adios2 import Adios, Stream
import matplotlib.pyplot as plt

def plot2d(Z, x_data, y_data, output_file='output.png', cmap='viridis'):
    """
    Plots a color plot of the data Z with x_data and y_data, saving it to a file.
    
    Parameters:
    - Z: 2D array, the data to plot
    - x_data: 1D array, the x-axis data
    - y_data: 1D array, the y-axis data
    - output_file: str, the filename to save the plot (default: 'output.png')
    - cmap: str, the colormap to use for the plot (default: 'viridis')
    """

    X, Y = np.meshgrid(x_data, y_data)
    x_label = "dim 1"
    y_label = "dim 0"
    title = "2d slices"
    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, Z, cmap=cmap)
    plt.colorbar(label='L or F value')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.savefig(output_file)
    plt.close()
    print(f"Image saved to {output_file}")


def read_data(varname, fr, start_coord, size_dims):
    data = fr.read(varname, start_coord, size_dims)
    # Read N-dim data but squeeze out the dimensions with size 1
    data = np.squeeze(data)
    return data


# Command-line arguments
if len(sys.argv) < 5:
    print("Usage: python3 slice2d.py BPfile  Variable  slice-index  x|y|z")
    sys.exit(1)

fname = sys.argv[1]
F_L = sys.argv[2]
slice_index = int(sys.argv[3])
slice_dim = sys.argv[4]

if slice_dim == 'x':
    var1 = 'y'
    var2 = 'z'
    start_coord = [slice_index, 0, 0]
elif slice_dim == 'y':
    var1 = 'x'
    var2 = 'z'
    start_coord = [0, slice_index, 0]
elif slice_dim == 'z':
    var1 = 'x'
    var2 = 'y'
    start_coord = [0, 0, slice_index]
else:
    print("slide_dim must be x/y/z")
    sys.exit(1)

# Reading data from file or stream as specified in adios2.xml
adios = Adios("adios2.xml")
io = adios.declare_io("WriteIO")
with Stream(io, fname, 'r') as f:
    for _ in f.steps():
        step = f.current_step()
        print(f"slice2d reads step {step}")
        var1_data = f.read(var1).flatten()
        var2_data = f.read(var2).flatten()
        vF = f.inquire_variable(F_L)
        shape = vF.shape()
        if slice_dim == 'x':
            size_dims = [1, shape[1], shape[2]]
        elif slice_dim == 'y':
            size_dims = [shape[0], 1, shape[2]]
        elif slice_dim == 'z':
            size_dims = [shape[0], shape[1], 1]
        data = f.read(F_L, start_coord, size_dims)

        # transform 3D array to 2D array (squeeze out the one with size 1)
        slice2D = np.squeeze(data)

        # Call the plot2d function and save the plot to a file
        plot2d(slice2D, var1_data, var2_data, output_file=f'{F_L}_{slice_dim}={slice_index}_{step}.png', cmap='plasma')