#!/bin/bash

qsub -A DD-20-28 -q qnvidia -l select=1:ncpus=16,walltime=4:00:00 -I
ml HDF5/1.10.4-gompi-2019.02-serial PGI CUDA/10.0.130
h5diff -v2 -p 0.000001 ../sampledata/sampleOutput.h5 step1Output.h5 /pos_x_final