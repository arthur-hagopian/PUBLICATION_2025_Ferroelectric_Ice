#!/usr/bin/env python3
# Arthur Hagopian <a.m.v.hagopian@lic.leidenuniv.nl>, version 05/01/2023


# THIS SCRIPT ALLOWS TO PLOT PROJECTED X-Y AVERAGES ON Z COORDINATE
# INPUT FILES ACCEPTED : chgcar.dat, locpot.dat, rhob.dat, rhoion.dat
# INPUT NEEDED : POSCAR


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from matplotlib.axes import Axes
from os import listdir, chdir
from os.path import isfile, isdir, join
import math
import sys
import matplotlib_parameters

# Get number of file(s) to plot
files_to_plot = input("Number of file(s) to plot : ")
files_to_plot = int(files_to_plot)

# Loop on the number of files / Get data
i = 0
files_names = []
list_z_values = []
list_projected_values = []
while i < files_to_plot:
    file_name = input("File number %s : " % (i+1))
    files_names.append(file_name)
    data = np.loadtxt(file_name)
    z_values = np.array(data[:, 0])
    projected_values = np.array(data[:, 1])
    proj_first = projected_values[0]
    projected_values = [proj - proj_first for proj in projected_values]
    list_z_values.append(z_values)
    list_projected_values.append(projected_values)
    i += 1

fig, ax = plt.subplots(1, figsize=(7, 5))

i = 0
while i < files_to_plot:
    plt.plot(list_z_values[i], list_projected_values[i], label=files_names[i])
    i += 1

plt.grid()
plt.legend()
plt.show()

