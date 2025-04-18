#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

# INPUT FILE NEEDED : "DOSCAR"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import sys

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

    return num_hydrogen, num_oxygen, hydrogen_positions, oxygen_positions

def read_bader(filename, num_hydrogen, num_oxygen):
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Extract charges
        charges = []
        i = 2 # Remove header
        while i < 2 + num_hydrogen + num_oxygen:
            charge = float(lines[i].split()[4])
            charges.append(charge)
            i += 1
    return charges

def compute_charges(charges, selected_hydrogen_groups, selected_oxygen_groups, labels, colors, alphas):
    fig, ax = plt.subplots()
    # Plot PDOS for selected atoms
    hydrogen_average_charges = []
    oxygen_average_charges = []
    for hydrogen_group in selected_hydrogen_groups:
        charge = 0
        hydrogen_charges = []
        oxygen_charges = []
        for hydrogen in hydrogen_group:
            charge += charges[hydrogen-1]
            hydrogen_charges.append(charges[hydrogen-1])
        charge = charge / len(hydrogen_group)
        hydrogen_average_charges.append(charge)
    for oxygen_group in selected_oxygen_groups:
        charge = 0
        for oxygen in oxygen_group:
            charge += charges[oxygen-1]
            oxygen_charges.append(charges[oxygen-1])
        charge = charge / len(oxygen_group)
        oxygen_average_charges.append(charge)
    print(hydrogen_average_charges)
    print(oxygen_average_charges)
    print(hydrogen_charges)
    print(oxygen_charges)

if __name__ == '__main__':
    # Specify the POSCAR & DOSCAR files
    poscar_file = 'CONTCAR'
    bader_file = 'ACF.dat'
    # Read the POSCAR file
    num_hydrogen, num_oxygen, hydrogen_positions, oxygen_positions = read_poscar(poscar_file)
    # Read the ACF.dat file
    charges = read_bader(bader_file, num_hydrogen, num_oxygen)
    # Define H atoms groups
    group_h_1 = []
    h_idx = 1
    for h in hydrogen_positions:
        if 0 < h[2] < 1.0:
            group_h_1.append(h_idx)
        h_idx += 1
    # Define O atoms groups
    group_o_1 = []
    o_idx = h_idx
    for o in oxygen_positions:
        if 0 < o[2] < 1.0:
            group_o_1.append(o_idx)
        o_idx += 1
    # Define groups
    selected_hydrogen_groups = [group_h_1]
    selected_oxygen_groups = [group_o_1]
    # Define for plot
    labels = ['Bilayer 1 : H', 'Bilayer 1 : O', 'Bilayer 2 : H', 'Bilayer 2 : O', 'Bilayer 3 : H', 'Bilayer 3 : O', 'Bilayer 4 : H', 'Bilayer 4 : O', 'Bilayer 5 : H', 'Bilayer 5 : O', 'Bilayer 6 : H', 'Bilayer 6 : O']
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'pink']
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    # Plot the summed PDOS
    compute_charges(charges, selected_hydrogen_groups, selected_oxygen_groups, labels, colors, alphas)

