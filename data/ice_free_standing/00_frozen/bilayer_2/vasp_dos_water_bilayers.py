#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

# INPUT FILE NEEDED : "DOSCAR"

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

    return hydrogen_positions, oxygen_positions

def parse_doscar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Extract header informations
        header = lines[:6]
        num_atoms = int(header[0].split()[0])
        efermi = float(header[5].split()[3])
        ndos = int(header[5].split()[2])
        # Extract total DOS
        total_dos_data = np.array([list(map(float,line.split())) for line in lines[6:6+ndos]])
        energies = total_dos_data[:,0] - efermi
        total_dos = total_dos_data[:,1]
        # Extract PDOS for each atom
        projected_dos = []
        start_idx = 6 + ndos + 1
        for atom in range(num_atoms):
            atom_dos_data = []
            for line in lines[start_idx:start_idx+ndos]:
                columns = list(map(float, line.split()))
                atom_dos_data.append(columns[1:])
            projected_dos.append(np.array(atom_dos_data))
            start_idx += ndos + 1

    return energies, total_dos, projected_dos, num_atoms

def plot_dos(energies, total_dos, projected_dos, selected_atoms_groups, labels, colors, alphas):
    fig, ax = plt.subplots()
    # Plot PDOS for selected atoms
    group = 0
    for selected_atoms in selected_atoms_groups:
        atom_dos_sum = []
        for atom_idx in selected_atoms:
            # Sum over all orbitals
            atom_dos = np.sum(projected_dos[atom_idx-1], axis=1)
            if len(atom_dos_sum) == 0:
                atom_dos_sum = atom_dos
            else:
                atom_dos_sum = atom_dos_sum + atom_dos
        atom_dos_sum = atom_dos_sum / len(selected_atoms)
        ax.plot(energies, atom_dos_sum, label=labels[group], color=colors[group], alpha=alphas[group])
        group += 1
    # Plot settings
    #ax.axvline(x=0, color='gray', linestyle='--', label='Fermi level')
    ax.set_xlabel("E - E$\mathrm{_F}$ (eV)")
    ax.set_ylabel("PDOS (a.u.)")
    ax.set_xlim([-15,5])
    ax.set_ylim(top=3.5)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax.legend()
    plt.show()

if __name__ == '__main__':
    # Specify the POSCAR & DOSCAR files
    poscar_file = 'POSCAR'
    doscar_file = 'DOSCAR'
    # Define limits between groups
    z_1 = 0.24
    # Read the POSCAR file
    hydrogen_positions, oxygen_positions = read_poscar(poscar_file)
    # Parse the DOSCAR file
    energies, total_dos, projected_dos, num_atoms = parse_doscar(doscar_file)
    # Define H atoms groups
    group_h_1 = []
    group_h_2 = []
    h_idx = 1
    for h in hydrogen_positions:
        if 0 < h[2] < z_1:
            group_h_1.append(h_idx)
        elif z_1 < h[2] < 1.0:
            group_h_2.append(h_idx)
        h_idx += 1
    # Define O atoms groups
    group_o_1 = []
    group_o_2 = []
    o_idx = h_idx
    for o in oxygen_positions:
        if 0 < o[2] < z_1:
            group_o_1.append(o_idx)
        elif z_1 < o[2] < 1.0:
            group_o_2.append(o_idx)
        o_idx += 1
    # Define groups
    selected_atoms_groups = [group_h_1, group_o_1, group_h_2, group_o_2]
    # Print to check
    print('H bilayer 1 : ' + str(selected_atoms_groups[0]))
    print('O bilayer 1 : ' + str(selected_atoms_groups[1]))
    print('H bilayer 2 : ' + str(selected_atoms_groups[2]))
    print('O bilayer 2 : ' + str(selected_atoms_groups[3]))
    # Define for plo
    labels = ['Bilayer 1 : H', 'Bilayer 1 : O', 'Bilayer 2 : H', 'Bilayer 2 : O']
    colors = ['blue', 'blue', 'red', 'red']
    alphas = [1.0, 0.5, 1.0, 0.5]
    # Plot the summed PDOS
    plot_dos(energies, total_dos, projected_dos, selected_atoms_groups, labels, colors, alphas)

