import sys
import numpy as np
import matplotlib.pyplot as plt
from adios2 import Adios, FileReader, Stream

if len(sys.argv) < 6:
    print(f"Usage: {sys.argv[0]} bpfile1 var1 bpfile2 var2 output_file")
    sys.exit(1)

fname1 = sys.argv[1]
var1 = sys.argv[2]
fname2 = sys.argv[3]
var2 = sys.argv[4]
output_file = sys.argv[5]

print(f"Reading data from the first file: {fname1}")
data1 = []
with Stream(fname1, "r") as f:
    for _ in f.steps():
        data1.append(f.read(var1))
print(f"Finished reading data from the first file: {fname1}")

print(f"Reading data from the second file: {fname2}")
data2 = []
with Stream(fname2, "r") as f:
    for _ in f.steps():
        data2.append(f.read(var2))
print(f"Finished reading data from the second file: {fname2}")

print("Converting lists to numpy arrays")
data1 = np.concatenate(data1)
data2 = np.concatenate(data2)
data1_flat = data1.flatten()
data2_flat = data2.flatten()


# match the sizes
min_length = min(len(data1_flat), len(data2_flat))
data1_flat = data1_flat[:min_length]
data2_flat = data2_flat[:min_length]

diff = data1_flat - data2_flat
print("Finished converting lists to numpy arrays")

print(f"Writing differences to the output file: {output_file}")
adios = Adios("adios2.xml")
ioOut = adios.declare_io("uncompressed error")
fout = Stream(ioOut, output_file, "w")

with Stream(output_file, "w") as s:
    fout.begin_step()
    fout.write("diff", diff, [len(diff)], [0], [len(diff)])
    fout.end_step()

fout.close()
print(f"Finished writing differences to the output file: {output_file}")

# Optionally, you can print out the differences or create a plot
print("Differences:")
print(diff)

# # If you want to create a plot
# plt.hist(diff, bins=50, alpha=0.5, label='Difference')
# plt.legend(loc='upper right')
# plt.title("Histogram of Differences")
# plt.xlabel("Difference")
# plt.ylabel("Frequency")
# plt.show()


