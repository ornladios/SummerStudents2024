import sys
import numpy as np
import matplotlib.pyplot as plt
from adios2 import Adios
from adios2 import FileReader, Stream

if len(sys.argv) < 8:
    print(f"Usage: {sys.argv[0]} bpfile funcname Lname bpfile1 funcname1 Lname1 output_file")
    sys.exit(1)

fname = sys.argv[1]
funcname = sys.argv[2]
Lname = sys.argv[3]
fname1 = sys.argv[4]
funcname1 = sys.argv[5]
Lname1 = sys.argv[6]
output_file = sys.argv[7]

print(f"Reading data from the first file: {fname}")
funcdata = []
Ldata = []
with Stream(fname, "r") as f:
    for _ in f.steps():
        funcdata.append(f.read(funcname))
        Ldata.append(f.read(Lname))
print(f"Finished reading data from the first file: {fname}")

print(f"Reading data from the second file: {fname1}")
funcdata1 = []
Ldata1 = []
with Stream(fname1, "r") as f:
    for _ in f.steps():
        funcdata1.append(f.read(funcname1))
        Ldata1.append(f.read(Lname1))
print(f"Finished reading data from the second file: {fname1}")

print("Converting lists to numpy arrays")
funcdata = np.concatenate(funcdata)
Ldata = np.concatenate(Ldata)
funcdata1 = np.concatenate(funcdata1)
Ldata1 = np.concatenate(Ldata1)
L = Ldata.flatten()
L1 = Ldata1.flatten()
F = funcdata.flatten()
F1 = funcdata1.flatten()

# Truncate the longer array to match the shorter one
min_length = min(len(L), len(L1))
L = L[:min_length]
L1 = L1[:min_length]

min_length = min(len(F), len(F1))
F = F[:min_length]
F1 = F1[:min_length]

Ldiff = L - L1
Fdiff = F - F1
print("Finished converting lists to numpy arrays")

print(f"Writing differences to the output file: {output_file}")
adios = Adios("adios2.xml")
ioOut = adios.declare_io("uncompressed error")
fout = Stream(ioOut, output_file, "w")

with Stream("output_file", "w") as s:
    fout.begin_step()
    fout.write("Fdiff", Fdiff, [len(Fdiff)], [0], [len(Fdiff)])
    fout.write("Ldiff", Ldiff, [len(Ldiff)], [0], [len(Ldiff)])
    fout.end_step()

fout.close()
print(f"Finished writing differences to the output file: {output_file}")

# Optionally, you can print out the differences or create a plot
print("Differences in F:")
print(Fdiff)
print("Differences in L:")
print(Ldiff)

# If you want to create a plot
# plt.hist(Fdiff, bins=50, alpha=0.5, label='Fdiff')
# plt.hist(Ldiff, bins=50, alpha=0.5, label='Ldiff')
# plt.legend(loc='upper right')
# plt.title("Histogram of Differences")
# plt.xlabel("Difference")
# plt.ylabel("Frequency")
# plt.show()

