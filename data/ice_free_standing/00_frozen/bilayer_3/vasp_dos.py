#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

# INPUT FILE NEEDED : "DOSCAR"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters

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

def plot_dos(energies, total_dos, projected_dos, selected_atoms):
    # Plot total DOS
    plt.subplot(211)
    plt.plot(energies, total_dos, label='Total DOS', color='black', linewidth='2')
    plt.ylabel("Total DOS (a.u.)")
    # Plot PDOS for selected atoms
    plt.subplot(212)
    for atom_idx in selected_atoms:
        # Sum over all orbitals
        atom_dos = np.sum(projected_dos[atom_idx-1], axis=1)
        plt.plot(energies, atom_dos, label=f'Atom {atom_idx}')
    # Plot settings
    #plt.axvline(x=0, color='gray', linestyle='--', label='Fermi level')
    plt.xlabel("E - E$\mathrm{_F}$ (eV)")
    plt.ylabel("PDOS (a.u.)")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    # Specify the DOSCAR file
    doscar_file = 'DOSCAR'
    # Parse the DOSCAR file
    energies, total_dos, projected_dos, num_atoms = parse_doscar(doscar_file)
    # Ask user for atom indices to plot PDOS
    selected_atoms = input(f'Enter the atom indices (1 to {num_atoms}) separated by spaces: ')
    selected_atoms = list(map(int,selected_atoms.split()))
    # Plot the DOS and PDOS
    plot_dos(energies, total_dos, projected_dos, selected_atoms)

