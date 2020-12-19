import h5py
import numpy as np

filename = "../sampledata/sampleOutput.h5"
filename1 = "step2Output.h5"

h5 = h5py.File(filename, 'r')
pos_x = h5['pos_x']  # VSTOXX futures data

h55 = h5py.File(filename1, 'r')
pos_xx = h55['pos_x_final']

print(pos_x[0])
print(pos_xx[0])

h5.close()