import sys
import numpy as np
import matplotlib.pyplot as plt
from adios2 import Stream

if len(sys.argv) < 6:
    print(f"Usage: {sys.argv[0]} need correct arguments")
    sys.exit(1)

fname = sys.argv[1]

varname = sys.argv[2]
fname1 = sys.argv[3]

varname1 = sys.argv[4]
step_ = int(sys.argv[5])  # Convert step_ to an integer

if int(np.log2(step_)) - np.log2(step_) != 0:
    print("Step size must be power of 2!!!!")
    sys.exit(1)

# Some empty lists to store the data

vardata = []

vardata1 = []

# Reading data from the first file
with Stream(fname, "r") as f:
    for _ in f.steps():
       
        vardata.append(f.read(varname))

# Reading data from the second file
with Stream(fname1, "r") as f:
    for _ in f.steps():
       
       vardata1.append(f.read(varname1))

dims_0 = vardata[0].shape

# Calculate the difference for only one timestep
for i in range(1):
    var = vardata[i][0:dims_0[0]:step_, 0:dims_0[1]:step_, 0:dims_0[2]:step_]
    var1 = vardata1[i][0:var.shape[0], 0:var.shape[1], 0:var.shape[2]]
    print(var.shape, var1.shape)
   
    assert len(var) == len(var1)

dims_1 = var.shape
vardiff = var - var1

# Visualize a single 2D slice from the 3D difference array
slice_idx = dims_1[0] // 2  # Middle slice for example
plt.figure()
plt.imshow(vardiff[slice_idx, :, :])
plt.colorbar()
plt.title(f"Slice {slice_idx} of Difference Array")
plt.savefig("difference_image.png")

rel_diff = np.abs(vardiff) / (var.max() - var.min())
print(np.where(np.reshape(rel_diff, [dims_1[0], dims_1[1], dims_1[2]]) > 0.1))

# Trim the ending points as needed
rel_diff2 = np.reshape(rel_diff, [dims_1[0], dims_1[1], dims_1[2]])[1:-2, 1:-2, 1:-2]
plt.figure()
plt.hist(rel_diff2.flatten(), bins=50)
plt.title("Histogram")
plt.show()
plt.savefig("difference_hist.png")
