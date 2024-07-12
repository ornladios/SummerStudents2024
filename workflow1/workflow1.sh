#!/bin/bash

EXE_DIR=..
SCRIPT_DIR=..

rm -rf *png *bp

mpirun -n 4 ${EXE_DIR}/calcF         F.bp 65 65 65 10 2           &
mpirun -n 2 ${EXE_DIR}/CalcLaplacian F.bp L.bp                    &
${EXE_DIR}/copier                    F.bp copyF.bp                &
python3 ${SCRIPT_DIR}/subtract.py    F.bp F x copyF.bp F diffF.bp   &
python3 ${SCRIPT_DIR}/slice.py       copyF.bp  F       32 x       &
python3 ${SCRIPT_DIR}/slice.py       L.bp      Laplace 32 x       &
python3 ${SCRIPT_DIR}/slice.py       diffF.bp  diff    32 x       &

wait 

echo "Completed workflow"
ls -ld *bp
ls -l *png
