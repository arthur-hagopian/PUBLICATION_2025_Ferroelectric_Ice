#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, updated 28/06/2025

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
import sys
import os

def parse_doscar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        header = lines[:6]
        num_atoms = int(header[0].split()[0])
        efermi = float(header[5].split()[3])
        ndos = int(header[5].split()[2])
        total_dos_data = np.array([list(map(float, line.split())) for line in lines[6:6+ndos]])
        energies = total_dos_data[:, 0] - efermi
        total_dos = total_dos_data[:, 1]
        mask = (energies <= 0)
        charge = np.trapz(total_dos[mask], energies[mask])
    return energies, total_dos, charge

def plot_multiple_dos(doscar_files):
    for i,file in enumerate(doscar_files):
        energies, total_dos, charge = parse_doscar(file)
        plt.plot(energies, total_dos, label=labels[i])
    
    plt.xlabel(r"$\mathrm{E - E_F\ (eV)}$")
    plt.ylabel("Total DOS (a.u.)")
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    doscar_files = ["1x1x1/DOSCAR", "2x2x1/DOSCAR", "3x3x1/DOSCAR", "4x4x1/DOSCAR", "5x5x1/DOSCAR"]
    labels = ["1x1x1", "2x2x1", "3x3x1", "4x4x1", "5x5x1"]

    for f in doscar_files:
        if not os.path.exists(f):
            print(f"File not found: {f}")
            sys.exit(1)

    plot_multiple_dos(doscar_files)

