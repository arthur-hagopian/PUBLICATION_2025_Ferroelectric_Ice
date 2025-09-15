#!/usr/bin/env python3
import os
import sys
import logging

import numpy as np
from scipy.interpolate import RegularGridInterpolator

# ——————— Logging setup ———————
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# === USER SETTINGS: just edit these three lines! ===
LOCPOT_FILE   = 'LOCPOT'                  # path to your LOCPOT file
INTERP_FACTOR = 1                         # integer >1 to up-sample grid
FRAC_POINT    = (0.5800782483774682,0.3846091147133421,0.1304954711713791)        # (x_frac, y_frac, z_frac)
# ===================================================

def read_locpot(filename):
    if not os.path.isfile(filename):
        logging.error(f"File {filename!r} not found.")
        sys.exit(1)
    with open(filename) as f:
        content = f.readlines()
    logging.info("LOCPOT file read successfully.")
    return content

def parse_lattice(content):
    scale = float(content[1].split()[0])
    lattice = np.array([[float(x) for x in content[2+i].split()] for i in range(3)])
    lattice *= scale
    logging.info(f"Parsed lattice (scale={scale}).")
    return lattice

def parse_grid_info(content):
    atom_counts = [int(x) for x in content[6].split()]
    atom_tot = sum(atom_counts)
    grid_idx = 9 + atom_tot
    NGX, NGY, NGZ = [int(x) for x in content[grid_idx].split()]
    logging.info(f"Grid dims: {NGX}×{NGY}×{NGZ}")
    return atom_tot, NGX, NGY, NGZ, grid_idx

def parse_potential(content, atom_tot, NGX, NGY, NGZ, grid_idx):
    total = NGX * NGY * NGZ
    vals = []
    i = grid_idx
    while len(vals) < total:
        i += 1
        vals.extend([float(x) for x in content[i].split()])
    if len(vals) != total:
        logging.error(f"Expected {total} vals, got {len(vals)}.")
        sys.exit(1)
    pot = np.zeros((NGX, NGY, NGZ), dtype=float)
    l = 0
    for k in range(NGZ):
        for j in range(NGY):
            for i in range(NGX):
                pot[i,j,k] = vals[l]
                l += 1
    logging.info("Potential grid parsed.")
    return pot

def interpolate_potential(pot, factor):
    if factor <= 1:
        return pot
    NGX, NGY, NGZ = pot.shape
    x = np.arange(NGX)/NGX
    y = np.arange(NGY)/NGY
    z = np.arange(NGZ)/NGZ
    interp = RegularGridInterpolator((x,y,z), pot,
                                     method='linear',
                                     bounds_error=False,
                                     fill_value=None)
    new_NX, new_NY, new_NZ = NGX*factor, NGY*factor, NGZ*factor
    fx = np.linspace(0,1,new_NX,endpoint=False)
    fy = np.linspace(0,1,new_NY,endpoint=False)
    fz = np.linspace(0,1,new_NZ,endpoint=False)
    ii, jj, kk = np.meshgrid(fx,fy,fz, indexing='ij')
    pts = np.vstack([ii.ravel(), jj.ravel(), kk.ravel()]).T
    vals = interp(pts)
    pot_fine = vals.reshape((new_NX,new_NY,new_NZ))
    logging.info(f"Interpolated to {new_NX}×{new_NY}×{new_NZ}")
    return pot_fine

def build_interpolator(pot):
    NGX, NGY, NGZ = pot.shape
    x = np.arange(NGX)/NGX
    y = np.arange(NGY)/NGY
    z = np.arange(NGZ)/NGZ
    return RegularGridInterpolator((x,y,z), pot,
                                   method='linear',
                                   bounds_error=False,
                                   fill_value=None)

def main():
    content = read_locpot(LOCPOT_FILE)
    lattice = parse_lattice(content)
    atom_tot, NGX, NGY, NGZ, grid_idx = parse_grid_info(content)
    pot = parse_potential(content, atom_tot, NGX, NGY, NGZ, grid_idx)

    if INTERP_FACTOR > 1:
        pot = interpolate_potential(pot, INTERP_FACTOR)
    logging.info(f"Using grid shape: {pot.shape}")

    interp_func = build_interpolator(pot)
    frac = np.array(FRAC_POINT, dtype=float) % 1.0

    # ——— FIXED LINE: convert to Python float ———
    value = interp_func(frac).item()
    print(f"Potential at fractional {tuple(frac)}: {value:.6f} eV")

if __name__ == '__main__':
    main()

