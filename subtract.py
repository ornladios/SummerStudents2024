import sys
import numpy as np
import matplotlib.pyplot as plt
from adios2 import FileReader, Stream


 

if len(sys.argv) < 7:
    print(f"Usage: {sys.argv[0]} bpfile funcname Lname bpfile1 funcname1 Lname1")
    sys.exit(1)

 

fname = sys.argv[1]
funcname = sys.argv[2]
Lname = sys.argv[3]
fname1 = sys.argv[4]
funcname1 = sys.argv[5]
Lname1 = sys.argv[6]

# some empty lists to store the data
funcdata = []
Ldata = []
funcdata1 = []
Ldata1 = []


# Reading data from the first file
with Stream(fname, "r") as f:
  for _ in f.steps():
    funcdata.append(f.read(funcname))
    Ldata.append(f.read(Lname))

# Reading data from the second file
with Stream(fname1, "r") as f:
    for _ in f.steps():
        funcdata1.append(f.read(funcname1))
        Ldata1.append(f.read(Lname1))


# Convert lists to numpy arrays
funcdata = np.concatenate(funcdata)
Ldata = np.concatenate(Ldata)
funcdata1 = np.concatenate(funcdata1)
Ldata1 = np.concatenate(Ldata1)
L = Ldata.flatten()
L1 = Ldata1.flatten()
F = funcdata.flatten()
F1 = funcdata1.flatten()
Ldiff = L - L1
Fdiff = F - F1

adios = Adios("adios2.xml")
ioOut = adios.declare_io("uncompressed error")
fout = Stream(ioOut, output_file, "w")

with Stream("output_file","w") as s:

        fout.begin_step()
        fout.write("Fdiff", Fdiff, [len(Fdiff)], [0], [len(Fdiff)])
        fout.write("Ldiff", Ldiff, [len(Ldiff)], [0], [len(Ldiff)])
        fout.end_step()

fout.close()

# for i in range(1):

#     N = len(Ldata[i])
#     n = len(Ldata1[i])
#     L = Ldata[i][0:N:2, 0:N:2, 0:N:2]
#     L1 = Ldata1[i][0:n]
#     print(L.shape, L1.shape)
#     # align? match
#     assert len(L) == len(L1)
#     # make sure the isze is 257 for big and 128 for small
#     # Ldiff = L -L1 fix
#     print(f"BIG {N} small {n}")


# L = L.flatten()
# L1 = L1.flatten()
# Ldiff = L - L1
# print(Ldiff[-1000:])
# #
# rel_diff = np.abs(Ldiff) / (L.max()-L.min())
# print(np.where(np.reshape(rel_diff, [129,129,129])>0.1))
# rel_diff2 = np.reshape(rel_diff, [129,129,129])[1:-2,1:-2,1:-2]
# plt.hist(rel_diff2.flatten())
# plt.title("Histogram ")
# plt.show()