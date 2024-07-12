import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from adios2 import Adios, Stream
import adios2.bindings as adios2

# Command-line arguments
if len(sys.argv) < 5:
    print("Usage: python3 histogram.py inputfile.bp var outputfilename.bp num_bin")
    sys.exit(1)

fname = sys.argv[1]
var = sys.argv[2]
output_name = sys.argv[3]
num_bin = int(sys.argv[4])



# Reading data from file or stream as specified in adios2.xml
adios = Adios("adios2.xml")
io = adios.declare_io("WriteIO")



with Stream(io, fname, 'r') as f:
    for _ in f.steps():
        step = f.current_step()
        print(f"reading step {step}")
        data = f.read(var)
        min_val = data.min()
        max_val = data.max()
        print(f"Variable {var} min = {min_val} max = {max_val}")
        
        # number of elements
        total_len = data.size 

        # Calculate bin size
        bin_size = (max_val - min_val) / num_bin

        # Initialize bins
        bins = np.zeros([num_bin])

        for value in data.flatten():
            # Determine the bin index
            bin_index = int((value - min_val) / bin_size)
            if bin_index == num_bin:  # Handle edge case where value is exactly max_val
                bin_index -= 1
            bins[bin_index] += 1

        # Print histogram
        # print(sum(bins), total_len)
        #
        # for i, val in enumerate(bins):
        #     print(f"Bin #{i} has {val} values.")

        # Validate the sum of bins equals the total length
        assert(sum(bins) == total_len)

        # fout.begin_step()
        # fout.write("bins", bins, [len(bins)], [0], [len(bins)])
        # fout.end_step()



# with Stream(io, output_name, 'r') as f:
#     for _ in f.steps():
#         bins_for_plot = f.read("bins")

# # Convert bins to numpy array
# bins_for_plot = np.array(bins_for_plot)

# Plotting the histogram
plt.bar(range(len(bins)), bins)
plt.xlabel('Bin Index')
plt.ylabel('Count')
plt.title('Histogram')
plt.show()
plt.savefig(f"{output_name}_histogram.png")


