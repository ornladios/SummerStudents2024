import sys
import numpy as np
from adios2 import Adios, Stream
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

FONTSIZE = 24


def plot2d(varname, Z, x_name, x_data, y_name, y_data, output_file='output.png', cmap='viridis'):
    """
    Plots a color plot of the data Z with x_data and y_data, saving it to a file.

    Parameters:
    - Z: 2D array, the data to plot
    - x_data: 1D array, the x-axis data
    - y_data: 1D array, the y-axis data
    - output_file: str, the filename to save the plot (default: 'output.png')
    - cmap: str, the colormap to use for the plot (default: 'viridis')
    """

    title = f"{x_name}-{y_name} slice of {varname}"

    gs = gridspec.GridSpec(1, 1)
    fig = plt.figure(1, figsize=(8, 8))
    ax = fig.add_subplot(gs[0, 0])
    minval = Z.min()
    maxval = Z.max()
    # print(f"Slice: Data array min = {minval}  max = {maxval}")
    colorax = ax.imshow(
        Z,
        origin="lower",
        interpolation="quadric",
        extent=[x_data.min(), x_data.max(), y_data.min(), y_data.max()],
        cmap=plt.get_cmap("gist_ncar"),
        vmin=minval, vmax=maxval
    )
    cbar = fig.colorbar(colorax, orientation="horizontal")
    cbar.ax.tick_params(labelsize=FONTSIZE-8)
    ax.set_title(title, fontsize=FONTSIZE)
    ax.set_xlabel(x_name, fontsize=FONTSIZE)
    ax.set_ylabel(y_name, fontsize=FONTSIZE)
    plt.tick_params(labelsize=FONTSIZE-8)

    plt.savefig(output_file)
    plt.close()
    print(f"Slice {varname}: Image saved to {output_file}")


def read_data(varname, fr, start_coord, size_dims):
    data = fr.read(varname, start_coord, size_dims)
    # Read N-dim data but squeeze out the dimensions with size 1
    data = np.squeeze(data)
    return data


# Command-line arguments
if len(sys.argv) < 5:
    print("Usage: python3 slice.py BPfile  Variable  slice-index  x|y|z")
    sys.exit(1)

fname = sys.argv[1]
varname = sys.argv[2]
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
    print("Slice: slide_dim must be x/y/z")
    sys.exit(1)

# Reading data from file or stream as specified in adios2.xml
adios = Adios("adios2.xml")
io = adios.declare_io("WriteIO")
with Stream(io, fname, 'r') as f:
    for _ in f.steps():
        step = f.current_step()
        print(f"Slice {fname}/{varname}: slice2d reads step {step}")

        v1 = f.inquire_variable(var1)
        v2 = f.inquire_variable(var2)
        vF = f.inquire_variable(varname)
        shape = vF.shape()
        if slice_dim == 'x':
            size_dims = [1, shape[1], shape[2]]
        elif slice_dim == 'y':
            size_dims = [shape[0], 1, shape[2]]
        elif slice_dim == 'z':
            size_dims = [shape[0], shape[1], 1]

        data = f.read(varname, start_coord, size_dims)

        # transform 3D array to 2D array (squeeze out the one with size 1)
        slice2D = np.squeeze(data)

        if v1:
            var1_data = f.read(var1).flatten()
        else:
            var1_data = np.arange(slice2D.shape[0])
            # print(f"data shape = {slice2D.shape}")
            # print(f"{var1} min = {var1_data.min()}  max = {var1_data.max()}")

        if v2:
            var2_data = f.read(var2).flatten()
        else:
            var2_data = np.arange(slice2D.shape[0])
            # print(f"{var2} min = {var2_data.min()}  max = {var2_data.max()}")

        # Call the plot2d function and save the plot to a file
        plot2d(varname, slice2D, var1, var1_data, var2, var2_data,
               output_file=f'{varname}_{slice_dim}={slice_index}_{step}.png',
               cmap='plasma')
