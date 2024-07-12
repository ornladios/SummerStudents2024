import sys
import numpy as np
from adios2 import Adios, Stream, bindings


if len(sys.argv) < 7:
    print(f"Usage: {sys.argv[0]} bpfile1 var1 varx bpfile2 var2 output_file")
    sys.exit(1)

fname1 = sys.argv[1]
var1 = sys.argv[2]
varx = sys.argv[3]
fname2 = sys.argv[4]
var2 = sys.argv[5]
output_file = sys.argv[6]

adios = Adios("adios2.xml")
io1 = adios.declare_io("WriteIO")
io2 = adios.declare_io("CopierOutput")
print(f"Subtract: open stream 1 {fname1}")
f1 = Stream(io1, fname1, "r")
print(f"Subtract: open stream 2 {fname2}")
f2 = Stream(io2, fname2, "r")
ioOut = adios.declare_io("uncompressed error")
print(f"Subtract: create output stream {output_file}")
fout = Stream(ioOut, output_file, "w")

step = 0
while True:
    print(f"Subtract: Step {step}")
    status = f1.begin_step()
    if status != bindings.StepStatus.OK:
        print(f"Subtract: No more steps or error reading first stream: {fname1}")
        break
    v1 = f1.inquire_variable(var1)
    vx = f1.inquire_variable(varx)
    shape1 = v1.shape()
    data1 = f1.read(v1)
    datax = f1.read(vx)
    step1 = f1.current_step()
    print(f"Subtract:    Data from stream 1, step = {step1} shape = {shape1}")
    f1.end_step()

    status = f2.begin_step()
    if status != bindings.StepStatus.OK:
        print(f"Subtract: No more steps or error reading second stream: {fname2}")
        break
    v2 = f2.inquire_variable(var2)
    shape2 = v2.shape()
    data2 = f2.read(v2)
    step2 = f2.current_step()
    print(f"Subtract:    Data from stream 2, step = {step2} shape = {shape2}")
    f2.end_step()

    if shape1 != shape2:
        print(f"Subtract: The shape of the two variables differ! {shape1}  and  {shape2}")
        break
    x = datax
    diff = data1 - data2
 
    fout.begin_step()
    start = np.zeros(3, dtype=np.int64)

    fout.write("diff", diff, [len(diff),len(diff),len(diff)], [0,0,0], [len(diff),len(diff),len(diff)])
    fout.write("x", x, [len(x)], [0], [len(x)])
    fout.end_step()

print("Subtract: Completed")
f1.close()
f2.close()
fout.close()

# # If you want to create a plot
# plt.hist(diff, bins=50, alpha=0.5, label='Difference')
# plt.legend(loc='upper right')
# plt.title("Histogram of Differences")
# plt.xlabel("Difference")
# plt.ylabel("Frequency")
# plt.show()


