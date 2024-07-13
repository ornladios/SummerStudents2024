# Workflow #1 

This workflow executes the following tasks:

- calcF  (producing F.bp, variable F)

three consumers of F:

- CalcLaplacioan (producing L.bp, variable Laplace)
- copier (producing copyF.bp, compressed with MGARD, variable F again)
- subtract (which also consumes copyF.bp), (producing diffF.bp, variable diff)

three 2D plotting:
- slice F
- slice Laplace
- slice diff

The streams are set to 
- use SST
- rendezvous set to block the producer until ALL consumers join
- block producer to make sure every consumer gets all steps


## Expected output:
Since the copier with compression is the slowest operation in this workflow, and all streams are set to be blocking, the order of the messages should be the same. 

```
$ ./workflow1.sh
Copier: Reading from: F.bp
Subtract: open stream 1 F.bp
Subtract: open stream 2 copyF.bp
CalcF: elapsed time = 0.00137377
Copier: Reading Step= 0
Copier: Writing Step= 0
Laplace: Max step: 10
Laplace: Delta T: 0.1
Laplace: output step: 0
Laplace: elapsed time = 0.000769377 step: 0
Slice L.bp/Laplace: slice2d reads step 0
Slice Laplace: Image saved to Laplace_x=32_0.png
Slice copyF.bp/F: slice2d reads step 0
Subtract: create output stream diffF.bp
Subtract: Step 0
Subtract:    Data from stream 1, step = 0 shape = [65, 65, 65]
CalcF: elapsed time = 3.31894
Copier: Reading Step= 1
Copier: Writing Step= 1
Laplace: output step: 1
Laplace: elapsed time = 0.00044775 step: 1
Slice L.bp/Laplace: slice2d reads step 1
Slice Laplace: Image saved to Laplace_x=32_1.png
Slice F: Image saved to F_x=32_0.png
Subtract:    Data from stream 2, step = 0 shape = [65, 65, 65]
Subtract: Step 0
Slice diffF.bp/diff: slice2d reads step 0
Subtract:    Data from stream 1, step = 1 shape = [65, 65, 65]
CalcF: elapsed time = 0.653612
Laplace: output step: 2
Laplace: elapsed time = 0.000405788 step: 2
Slice L.bp/Laplace: slice2d reads step 2
Slice Laplace: Image saved to Laplace_x=32_2.png
Slice diff: Image saved to diff_x=32_0.png
Copier: Reading Step= 2
Slice copyF.bp/F: slice2d reads step 1
Copier: Writing Step= 2
Subtract:    Data from stream 2, step = 1 shape = [65, 65, 65]
Subtract: Step 0
Slice diffF.bp/diff: slice2d reads step 1
Subtract:    Data from stream 1, step = 2 shape = [65, 65, 65]
CalcF: elapsed time = 2.89775
Laplace: output step: 3
Laplace: elapsed time = 0.000432968 step: 3
Slice L.bp/Laplace: slice2d reads step 3
Slice F: Image saved to F_x=32_1.png
Slice Laplace: Image saved to Laplace_x=32_3.png
Slice diff: Image saved to diff_x=32_1.png
Slice copyF.bp/F: slice2d reads step 2
Copier: Reading Step= 3
Copier: Writing Step= 3
Subtract:    Data from stream 2, step = 2 shape = [65, 65, 65]
Subtract: Step 0
Slice diffF.bp/diff: slice2d reads step 2
Subtract:    Data from stream 1, step = 3 shape = [65, 65, 65]
CalcF: elapsed time = 3.09329
Laplace: output step: 4
Laplace: elapsed time = 0.000595808 step: 4
Slice L.bp/Laplace: slice2d reads step 4
Slice F: Image saved to F_x=32_2.png
Slice Laplace: Image saved to Laplace_x=32_4.png
Slice diff: Image saved to diff_x=32_2.png
Copier: Reading Step= 4
Slice copyF.bp/F: slice2d reads step 3
Copier: Writing Step= 4
Subtract:    Data from stream 2, step = 3 shape = [65, 65, 65]
Subtract: Step 0
Subtract:    Data from stream 1, step = 4 shape = [65, 65, 65]
CalcF: elapsed time = 3.33332
Laplace: output step: 5
Laplace: elapsed time = 0.000864983 step: 5
Slice L.bp/Laplace: slice2d reads step 5
Slice diffF.bp/diff: slice2d reads step 3
Slice Laplace: Image saved to Laplace_x=32_5.png
Slice diff: Image saved to diff_x=32_3.png
Slice F: Image saved to F_x=32_3.png
Slice copyF.bp/F: slice2d reads step 4
Copier: Reading Step= 5
Copier: Writing Step= 5
Subtract:    Data from stream 2, step = 4 shape = [65, 65, 65]
Subtract: Step 0
Slice diffF.bp/diff: slice2d reads step 4
Subtract:    Data from stream 1, step = 5 shape = [65, 65, 65]
CalcF: elapsed time = 3.05105
Laplace: output step: 6
Laplace: elapsed time = 0.000554323 step: 6
Slice L.bp/Laplace: slice2d reads step 6
Slice F: Image saved to F_x=32_4.png
Slice Laplace: Image saved to Laplace_x=32_6.png
Slice diff: Image saved to diff_x=32_4.png
Slice copyF.bp/F: slice2d reads step 5
Copier: Reading Step= 6
Copier: Writing Step= 6
Subtract:    Data from stream 2, step = 5 shape = [65, 65, 65]
Subtract: Step 0
Slice diffF.bp/diff: slice2d reads step 5
Subtract:    Data from stream 1, step = 6 shape = [65, 65, 65]
CalcF: elapsed time = 3.22865
Laplace: output step: 7
Laplace: elapsed time = 0.000804424 step: 7
Slice L.bp/Laplace: slice2d reads step 7
Slice F: Image saved to F_x=32_5.png
Slice Laplace: Image saved to Laplace_x=32_7.png
Slice diff: Image saved to diff_x=32_5.png
Slice copyF.bp/F: slice2d reads step 6
Copier: Reading Step= 7
Copier: Writing Step= 7
Subtract:    Data from stream 2, step = 6 shape = [65, 65, 65]
Subtract: Step 0
Slice diffF.bp/diff: slice2d reads step 6
Subtract:    Data from stream 1, step = 7 shape = [65, 65, 65]
CalcF: elapsed time = 3.21286
Laplace: output step: 8
Laplace: elapsed time = 0.000651836 step: 8
Slice L.bp/Laplace: slice2d reads step 8
Slice F: Image saved to F_x=32_6.png
Slice diff: Image saved to diff_x=32_6.png
Slice Laplace: Image saved to Laplace_x=32_8.png
Copier: Reading Step= 8
Slice copyF.bp/F: slice2d reads step 7
Copier: Writing Step= 8
Subtract:    Data from stream 2, step = 7 shape = [65, 65, 65]
Subtract: Step 0
Slice diffF.bp/diff: slice2d reads step 7
Subtract:    Data from stream 1, step = 8 shape = [65, 65, 65]
CalcF: elapsed time = 3.16804
Laplace: output step: 9
Laplace: elapsed time = 0.000432253 step: 9
Slice L.bp/Laplace: slice2d reads step 9
Slice diff: Image saved to diff_x=32_7.png
Slice F: Image saved to F_x=32_7.png
Slice Laplace: Image saved to Laplace_x=32_9.png
Slice copyF.bp/F: slice2d reads step 8
Copier: Reading Step= 9
Copier: Writing Step= 9
Subtract:    Data from stream 2, step = 8 shape = [65, 65, 65]
Subtract: Step 0
Slice diffF.bp/diff: slice2d reads step 8
Subtract:    Data from stream 1, step = 9 shape = [65, 65, 65]
CalcF: Completed. Output is F.bp
Slice F: Image saved to F_x=32_8.png
Slice diff: Image saved to diff_x=32_8.png
Slice copyF.bp/F: slice2d reads step 9
Subtract:    Data from stream 2, step = 9 shape = [65, 65, 65]
Subtract: Step 0
Subtract: No more steps or error reading first stream: F.bp
Subtract: Completed
Slice diffF.bp/diff: slice2d reads step 9
Slice F: Image saved to F_x=32_9.png
Slice diff: Image saved to diff_x=32_9.png
Completed workflow
ls: cannot access '*bp': No such file or directory
-rw-r--r-- 1 adios adios  28353 Jul 12 15:23 'F_x=32_0.png'
-rw-r--r-- 1 adios adios  41673 Jul 12 15:23 'F_x=32_1.png'
-rw-r--r-- 1 adios adios  53411 Jul 12 15:24 'F_x=32_2.png'
-rw-r--r-- 1 adios adios  58934 Jul 12 15:24 'F_x=32_3.png'
-rw-r--r-- 1 adios adios  65530 Jul 12 15:24 'F_x=32_4.png'
-rw-r--r-- 1 adios adios  67497 Jul 12 15:24 'F_x=32_5.png'
-rw-r--r-- 1 adios adios  68322 Jul 12 15:24 'F_x=32_6.png'
-rw-r--r-- 1 adios adios  69334 Jul 12 15:24 'F_x=32_7.png'
-rw-r--r-- 1 adios adios  69595 Jul 12 15:24 'F_x=32_8.png'
-rw-r--r-- 1 adios adios  70314 Jul 12 15:24 'F_x=32_9.png'
-rw-r--r-- 1 adios adios  58518 Jul 12 15:23 'Laplace_x=32_0.png'
-rw-r--r-- 1 adios adios  61381 Jul 12 15:23 'Laplace_x=32_1.png'
-rw-r--r-- 1 adios adios  67077 Jul 12 15:23 'Laplace_x=32_2.png'
-rw-r--r-- 1 adios adios  69562 Jul 12 15:23 'Laplace_x=32_3.png'
-rw-r--r-- 1 adios adios  72958 Jul 12 15:24 'Laplace_x=32_4.png'
-rw-r--r-- 1 adios adios  75907 Jul 12 15:24 'Laplace_x=32_5.png'
-rw-r--r-- 1 adios adios  79587 Jul 12 15:24 'Laplace_x=32_6.png'
-rw-r--r-- 1 adios adios  83042 Jul 12 15:24 'Laplace_x=32_7.png'
-rw-r--r-- 1 adios adios  84866 Jul 12 15:24 'Laplace_x=32_8.png'
-rw-r--r-- 1 adios adios  86496 Jul 12 15:24 'Laplace_x=32_9.png'
-rw-r--r-- 1 adios adios 197650 Jul 12 15:23 'diff_x=32_0.png'
-rw-r--r-- 1 adios adios 197933 Jul 12 15:23 'diff_x=32_1.png'
-rw-r--r-- 1 adios adios 186360 Jul 12 15:24 'diff_x=32_2.png'
-rw-r--r-- 1 adios adios 187199 Jul 12 15:24 'diff_x=32_3.png'
-rw-r--r-- 1 adios adios 190065 Jul 12 15:24 'diff_x=32_4.png'
-rw-r--r-- 1 adios adios 199271 Jul 12 15:24 'diff_x=32_5.png'
-rw-r--r-- 1 adios adios 203561 Jul 12 15:24 'diff_x=32_6.png'
-rw-r--r-- 1 adios adios 175686 Jul 12 15:24 'diff_x=32_7.png'
-rw-r--r-- 1 adios adios 189853 Jul 12 15:24 'diff_x=32_8.png'
-rw-r--r-- 1 adios adios 170047 Jul 12 15:24 'diff_x=32_9.png'
```
