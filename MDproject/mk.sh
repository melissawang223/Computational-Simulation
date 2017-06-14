#!/bin/sh

source /usr/usc/openmpi/1.6.4/setup.sh

rm -vf gr pmd
gcc -lm gr.c -o gr
mpicc -w pmd.c -o pmd
