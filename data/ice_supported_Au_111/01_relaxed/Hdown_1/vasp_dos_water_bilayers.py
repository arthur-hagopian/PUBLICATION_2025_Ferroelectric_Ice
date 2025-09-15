#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

# INPUT FILE NEEDED : "DOSCAR"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def read_poscar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Extract header informations
        header = lines[:9]
        num_gold = int(header[6].split()[0])
        num_hydrogen = int(header[6].split()[1])
        num_oxygen = int(header[6].split()[2])
        # Extract positions
        gold_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9:9+num_gold]])
        hydrogen_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9+num_gold:9+num_gold+num_hydrogen]])
        oxygen_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9+num_gold+num_hydrogen:9+num_gold+num_hydrogen+num_oxygen]])

    return gold_positions, hydrogen_positions, oxygen_positions

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

def plot_dos(energies, total_dos, projected_dos, selected_atoms_groups, labels, colors):
    fig, ax = plt.subplots()
    group = 0
    for index, selected_atoms in enumerate(selected_atoms_groups):
        if not selected_atoms:  # Skip if the list is empty
            group += 1
            continue

        atom_dos_sum = None
        for atom_idx in selected_atoms:
            atom_dos = np.sum(projected_dos[atom_idx - 1], axis=1)
            if atom_dos_sum is None:
                atom_dos_sum = atom_dos
            else:
                atom_dos_sum += atom_dos

        atom_dos_sum /= len(selected_atoms)

        linestyle = '--' if index % 2 == 0 else '-'
        alpha = 0.7 if index % 2 == 0 else 0.9
        linewidth = 2.0 if index % 2 == 0 else 3.0

        ax.plot(
            energies,
            atom_dos_sum,
            label=labels[group] if index % 2 != 0 else None,  # Only label O for legend
            color=colors[group],
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha
        )
        group += 1

    # Plot settings
    ax.set_xlabel("E - E$\mathrm{_F}$ (eV)")
    ax.set_ylabel("PDOS (a.u.)")
    ax.set_xlim([-15, 5])
    ax.set_ylim(top=2.5)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    #ax.legend(ncol=6)

    # Save figure
    output_folder = "/home/ahagopian/WORK/postdoc_leiden/projects/06_Ice_XI/images/"
    current_path = os.getcwd()
    parts = current_path.split(os.sep)[-3:]  # take last 3 parts
    fig_name = "_".join(parts)  # join with underscore because / is illegal in filenames
    output_path = os.path.join(output_folder, f"{fig_name}.svg")
    plt.savefig(output_path, format='svg', bbox_inches='tight')

    plt.show()

if __name__ == '__main__':
    # Specify the POSCAR & DOSCAR files
    poscar_file = 'CONTCAR'
    doscar_file = 'DOSCAR'
    # Define limits between groups
    z_1 = 0.24
    z_2 = 0.32
    z_3 = 0.40
    z_4 = 0.48
    z_5 = 0.56
    # Read the POSCAR file
    gold_positions, hydrogen_positions, oxygen_positions = read_poscar(poscar_file)
    # Parse the DOSCAR file
    energies, total_dos, projected_dos, num_atoms = parse_doscar(doscar_file)
    # Skip Au group
    au_idx = len(gold_positions)
    # Define H atoms groups
    group_h_1 = []
    group_h_2 = []
    group_h_3 = []
    group_h_4 = []
    group_h_5 = []
    group_h_6 = []
    h_idx = au_idx + 1
    for h in hydrogen_positions:
        if 0 < h[2] < z_1:
            group_h_1.append(h_idx)
        elif z_1 < h[2] < z_2:
            group_h_2.append(h_idx)
        elif z_2 < h[2] < z_3:
            group_h_3.append(h_idx)
        elif z_3 < h[2] < z_4:
            group_h_4.append(h_idx)
        elif z_4 < h[2] < z_5:
            group_h_5.append(h_idx)
        elif z_5 < h[2] < 1.0:
            group_h_6.append(h_idx)
        h_idx += 1
    # Define O atoms groups
    group_o_1 = []
    group_o_2 = []
    group_o_3 = []
    group_o_4 = []
    group_o_5 = []
    group_o_6 = []
    o_idx = h_idx
    for o in oxygen_positions:
        if 0 < o[2] < z_1:
            group_o_1.append(o_idx)
        elif z_1 < o[2] < z_2:
            group_o_2.append(o_idx)
        elif z_2 < o[2] < z_3:
            group_o_3.append(o_idx)
        elif z_3 < o[2] < z_4:
            group_o_4.append(o_idx)
        elif z_4 < o[2] < z_5:
            group_o_5.append(o_idx)
        elif z_5 < o[2] < 1.0:
            group_o_6.append(o_idx)
        o_idx += 1
    # Define groups
    selected_atoms_groups = [group_h_1, group_o_1, group_h_2, group_o_2, group_h_3, group_o_3, group_h_4, group_o_4, group_h_5, group_o_5, group_h_6, group_o_6]
    # Define for plot
    labels = ['Bilayer 1', 'Bilayer 1', 'Bilayer 2', 'Bilayer 2', 'Bilayer 3', 'Bilayer 3', 'Bilayer 4', 'Bilayer 4', 'Bilayer 5', 'Bilayer 5', 'Bilayer 6', 'Bilayer 6']
    colors = [
    "#377eb8", "#377eb8",  # Bilayer 1 H, O
    "#e41a1c", "#e41a1c",  # Bilayer 2 H, O
    "#4daf4a", "#4daf4a",  # Bilayer 3 H, O
    "#ff7f00", "#ff7f00",  # Bilayer 4 H, O
    "#984ea3", "#984ea3",  # Bilayer 5 H, O
    "#a65628", "#a65628"   # Bilayer 6 H, O
    ]
    # Plot the summed PDOS
    plot_dos(energies, total_dos, projected_dos, selected_atoms_groups, labels, colors)
