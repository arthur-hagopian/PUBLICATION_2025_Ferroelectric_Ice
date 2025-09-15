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

fig, ax = plt.subplots()

files_to_plot = ["locpot_1.dat",
                 "locpot_2.dat",
                 "locpot_3.dat",
                 "locpot_4.dat",
                 "locpot_5.dat",
                 "locpot_6.dat",
                 ]

colors = [
        "#377eb8",
        "#e41a1c",
        "#4daf4a",
        "#ff7f00",
        "#984ea3",
        "#a65628",
        ]

labels = ['Bilayer 1','Bilayer 2','Bilayer 3','Bilayer 4','Bilayer 5','Bilayer 6']

# LOOP
for index,file in enumerate(files_to_plot):
    data = np.loadtxt(file)
    z_values = np.array(data[:, 0])
    projected_values = np.array(data[:, 1])
    proj_first = projected_values[0]
    projected_values = [proj - proj_first for proj in projected_values]
    if index == 0 or index == 1:
        plt.plot(z_values, projected_values, label=labels[index], color=colors[index])
    else:
        plt.plot(z_values, projected_values, label=labels[index], color=colors[index], linewidth=1.5, linestyle='--')


plt.xlim(4,33)
plt.xlabel("Z-axis (Å)")
plt.ylabel("Electrostatic potential (eV)")
plt.grid(False)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3)
plt.show()

