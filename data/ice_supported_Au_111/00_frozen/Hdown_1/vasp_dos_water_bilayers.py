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
        num_gold = int(header[6].split()[0])
        num_hydrogen = int(header[6].split()[1])
        num_oxygen = int(header[6].split()[2])
        # Extract positions
        gold_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9:9+num_gold]])
        hydrogen_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9+num_gold:9+num_gold+num_hydrogen]])
        oxygen_positions = np.array([list(map(float,line.split()[:3])) for line in lines[9+num_gold+num_hydrogen:9+num_gold+num_hydrogen+num_oxygen]])

    return gold_positions,hydrogen_positions, oxygen_positions

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
            k = []
            for line in lines[start_idx:start_idx+ndos]:
                columns = list(map(float, line.split()))
                atom_dos_data.append(columns[1:])
                k.append(sum(columns[1:]))
            projected_dos.append(np.array(atom_dos_data))
            mask = (energies <= 0)
            energies_for_int = energies[mask]
            k = np.array(k)
            atom_dos_sum_for_int = k[mask]
            charge = np.trapz(energies_for_int,atom_dos_sum_for_int)
            print(charge)
            start_idx += ndos + 1

    return energies, total_dos, projected_dos, num_atoms

def plot_dos(energies, total_dos, projected_dos, selected_atoms_groups, labels, colors, alphas):
    fig, ax = plt.subplots()
    # Plot PDOS for selected atoms
    group = 0
    #ax.plot(energies, total_dos, label='Total DOS', color='black', linewidth='2')
    atom_dos_sum_all = []
    for selected_atoms in selected_atoms_groups:
        atom_dos_sum = []
        for atom_idx in selected_atoms:
            # Sum over all orbitals
            atom_dos = np.sum(projected_dos[atom_idx-1], axis=1)
            if len(atom_dos_sum) == 0:
                atom_dos_sum = atom_dos
            else:
                atom_dos_sum = atom_dos_sum + atom_dos
        atom_dos_sum_all.append(atom_dos_sum)
        atom_dos_sum = atom_dos_sum / len(selected_atoms)
        # Integral total DOS till efermi
        mask = (energies <= 0)
        energies_for_int = energies[mask]
        atom_dos_sum_for_int = atom_dos_sum[mask]
        charge = np.trapz(energies_for_int,atom_dos_sum_for_int)
        print(f'Total charge {labels[group]} : {charge}')
        # Plot
        ax.plot(energies, atom_dos_sum, label=labels[group], color=colors[group], alpha=alphas[group], linestyle='-')
        group += 1
    sum_dos = atom_dos_sum_all[0] + atom_dos_sum_all[1]
    #ax.plot(energies, sum_dos, label='$\mathrm{\sum}$ PDOS', color=colors[0], alpha=1.0)
    # Plot settings
    #ax.axvline(x=0, color='gray', linestyle='--', label='Fermi level')
    ax.set_xlabel("E - E$\mathrm{_F}$ (eV)")
    ax.set_ylabel("PDOS (a.u.)")
    ax.set_xlim([-15,5])
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax.legend()
    plt.show()

if __name__ == '__main__':
    # Specify the POSCAR & DOSCAR files
    poscar_file = 'POSCAR'
    doscar_file = 'DOSCAR'
    # Read the POSCAR file
    gold_positions, hydrogen_positions, oxygen_positions = read_poscar(poscar_file)
    # Parse the DOSCAR file
    energies, total_dos, projected_dos, num_atoms = parse_doscar(doscar_file)
    # Skip Au group
    au_idx = len(gold_positions)
    # Define H atoms groups
    group_h_1 = []
    h_idx = au_idx + 1
    for h in hydrogen_positions:
        group_h_1.append(h_idx)
        h_idx += 1
    # Define O atoms groups
    group_o_1 = []
    o_idx = h_idx
    for o in oxygen_positions:
        group_o_1.append(o_idx)
        o_idx += 1
    # Define groups
    selected_atoms_groups = [group_h_1, group_o_1]
    # Print to check
    print('H bilayer 1 : ' + str(selected_atoms_groups[0]))
    print('O bilayer 1 : ' + str(selected_atoms_groups[1]))
    # Define for plot
    labels = ['Bilayer 1 : H', 'Bilayer 1 : O']
    colors = ['blue', 'blue']
    alphas = [1.0, 0.5]
    # Plot the summed PDOS
    plot_dos(energies, total_dos, projected_dos, selected_atoms_groups, labels, colors, alphas)

