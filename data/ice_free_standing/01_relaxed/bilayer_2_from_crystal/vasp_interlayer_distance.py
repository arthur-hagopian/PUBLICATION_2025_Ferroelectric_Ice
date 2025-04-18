#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

# INPUT FILE NEEDED : "CONTCAR"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def read_poscar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Extract header informations
        header = lines[:9]
        num_hydrogen = int(header[6].split()[0])
        num_oxygen = int(header[6].split()[1])
        # Extract positions
        hydrogen_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9:9+num_hydrogen]])
        oxygen_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9+num_hydrogen:9+num_hydrogen+num_oxygen]])
        all_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9:9+num_hydrogen+num_oxygen]])

    return hydrogen_positions, oxygen_positions, all_positions

if __name__ == '__main__':
    # Specify the POSCAR & DOSCAR files
    poscar_file = 'CONTCAR'
    # Define limits between groups
    z_1_2 = 0.24
    # Read the POSCAR file
    hydrogen_positions, oxygen_positions, all_positions = read_poscar(poscar_file)
    # Define H atoms groups
    group_1 = []
    group_2 = []
    for x in all_positions:
        if 0 < x[2] < z_1_2:
            group_1.append(x[2])
        elif z_1_2 < x[2] < 1.0:
            group_2.append(x[2])
    # Average
    z_group_1 = sum(group_1)/len(group_1)
    z_group_2 = sum(group_2)/len(group_2)
    print(z_group_1,z_group_2)
    # Compute interlayer distance
    int_dist_1_2 = z_group_2 - z_group_1
    print(f'Interlayer distance : {int_dist_1_2}')

