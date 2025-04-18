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
    z_2_3 = 0.32
    z_3_4 = 0.40
    z_4_5 = 0.48
    z_5_6 = 0.56
    # Read the POSCAR file
    hydrogen_positions, oxygen_positions, all_positions = read_poscar(poscar_file)
    # Define H atoms groups
    group_1 = []
    group_2 = []
    group_3 = []
    group_4 = []
    group_5 = []
    group_6 = []
    for x in all_positions:
        if 0 < x[2] < z_1_2:
            group_1.append(x[2])
        elif z_1_2 < x[2] < z_2_3:
            group_2.append(x[2])
        elif z_2_3 < x[2] < z_3_4:
            group_3.append(x[2])
        elif z_3_4 < x[2] < z_4_5:
            group_4.append(x[2])
        elif z_4_5 < x[2] < z_5_6:
            group_5.append(x[2])
        elif z_5_6 < x[2] < 1.0:
            group_6.append(x[2])
    # Average
    z_group_1 = sum(group_1)/len(group_1)
    z_group_2 = sum(group_2)/len(group_2)
    z_group_3 = sum(group_3)/len(group_3)
    z_group_4 = sum(group_4)/len(group_4)
    z_group_5 = sum(group_5)/len(group_5)
    z_group_6 = sum(group_6)/len(group_6)
    print(z_group_1,z_group_2,z_group_3,z_group_4,z_group_5,z_group_6)
    # Compute interlayer distance
    int_dist_1_2 = z_group_2 - z_group_1
    int_dist_2_3 = z_group_3 - z_group_2
    int_dist_3_4 = z_group_4 - z_group_3
    int_dist_4_5 = z_group_5 - z_group_4
    int_dist_5_6 = z_group_6 - z_group_5
    print(f'Interlayer distance : {int_dist_1_2}')
    print(f'Interlayer distance : {int_dist_2_3}')
    print(f'Interlayer distance : {int_dist_3_4}')
    print(f'Interlayer distance : {int_dist_4_5}')
    print(f'Interlayer distance : {int_dist_5_6}')
