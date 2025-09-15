#!/usr/bin/env python3
import os
import sys
import logging
import argparse

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ——————— Logging setup ———————
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def read_locpot(filename):
    if not os.path.isfile(filename):
        logging.error(f"File {filename} not found.")
        sys.exit(1)
    with open(filename) as f:
        content = f.readlines()
    logging.info("LOCPOT file read successfully.")
    return content


def parse_lattice(content):
    scale = float(content[1].split()[0])
    lat = np.array([[float(x) for x in content[2+i].split()] for i in range(3)])
    lat *= scale
    logging.info(f"Parsed lattice (scale={scale}):\n{lat}")
    return lat


def parse_grid_info(content):
    atom_counts = [int(x) for x in content[6].split()]
    atom_tot = sum(atom_counts)
    grid_idx = 9 + atom_tot
    NGX, NGY, NGZ = [int(x) for x in content[grid_idx].split()]
    logging.info(f"Grid dims: NGX={NGX}, NGY={NGY}, NGZ={NGZ}")
    return atom_tot, NGX, NGY, NGZ, grid_idx


def parse_potential(content, atom_tot, NGX, NGY, NGZ, grid_idx):
    total = NGX * NGY * NGZ
    vals = []
    i = grid_idx
    while len(vals) < total:
        i += 1
        vals.extend([float(x) for x in content[i].split()])
    if len(vals) != total:
        logging.error(f"Expected {total} values, got {len(vals)}")
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
    NGX, NGY, NGZ = pot.shape
    x = np.arange(NGX)/NGX
    y = np.arange(NGY)/NGY
    z = np.arange(NGZ)/NGZ
    interp = RegularGridInterpolator((x,y,z), pot,
                                     method='linear',
                                     bounds_error=False,
                                     fill_value=None)
    new_NX, new_NY, new_NZ = NGX*factor, NGY*factor, NGZ*factor
    fx = np.linspace(0, 1, new_NX, endpoint=False)
    fy = np.linspace(0, 1, new_NY, endpoint=False)
    fz = np.linspace(0, 1, new_NZ, endpoint=False)
    ii, jj, kk = np.meshgrid(fx, fy, fz, indexing='ij')
    pts = np.vstack([ii.ravel(), jj.ravel(), kk.ravel()]).T
    vals = interp(pts)
    pot_fine = vals.reshape((new_NX, new_NY, new_NZ))
    logging.info(f"Interpolated grid to ({new_NX}, {new_NY}, {new_NZ})")
    return pot_fine


def make_cartesian_coords(lattice, NGX, NGY, NGZ):
    i = np.arange(NGX)
    j = np.arange(NGY)
    k = np.arange(NGZ)
    ii, jj, kk = np.meshgrid(i, j, k, indexing='ij')
    frac = np.vstack([
        ii.ravel()/NGX,
        jj.ravel()/NGY,
        kk.ravel()/NGZ
    ]).T
    return frac @ lattice

def plot_3d_potential(lattice, pot, sample=4, vmin=None, vmax=None):
    NGX, NGY, NGZ = pot.shape
    coords = make_cartesian_coords(lattice, NGX, NGY, NGZ)
    vals   = pot.ravel()

    # 1) down‐sample mask
    idx = np.indices((NGX,NGY,NGZ))
    samp_mask = (
        (idx[0].ravel() % sample == 0) &
        (idx[1].ravel() % sample == 0) &
        (idx[2].ravel() % sample == 0)
    )
    c_ds = coords[samp_mask]
    v_ds = vals[samp_mask]

    # 2) mask out‐of‐range values
    if vmin is not None:
        in_range = v_ds >= vmin
    else:
        in_range = np.ones_like(v_ds, dtype=bool)
    if vmax is not None:
        in_range &= v_ds <= vmax

    c_plot = c_ds[in_range]
    v_plot = v_ds[in_range]

    # 3) scatter only the in‐range points
    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111, projection='3d')
    sc  = ax.scatter(
        c_plot[:,0], c_plot[:,1], c_plot[:,2],
        c=v_plot, cmap='seismic',
        vmin=vmin, vmax=vmax,
        s=3, alpha=0.8
    )
    cbar = fig.colorbar(sc, ax=ax, label='Potential (eV)')

    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('3D LOCPOT (masked outside range)')

    # 4) equal‐aspect
    spans = c_plot.max(axis=0) - c_plot.min(axis=0)
    mr = spans.max()/2
    mid = (c_plot.max(axis=0) + c_plot.min(axis=0)) / 2
    ax.set_xlim(mid[0]-mr, mid[0]+mr)
    ax.set_ylim(mid[1]-mr, mid[1]+mr)
    ax.set_zlim(mid[2]-mr, mid[2]+mr)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot VASP LOCPOT in 3D with optional interpolation and value range."
    )
    parser.add_argument('file', help='Path to LOCPOT file')
    parser.add_argument('--interp', '-i', type=int, default=1,
                        help='Interpolation factor (integer)')
    parser.add_argument('--sample', '-s', type=int, default=3,
                        help='Sampling rate for scatter plot')
    parser.add_argument('--vmin', type=float, default=None,
                        help='Minimum potential for colormap')
    parser.add_argument('--vmax', type=float, default=None,
                        help='Maximum potential for colormap')
    args = parser.parse_args()

    content = read_locpot(args.file)
    lattice = parse_lattice(content)
    atom_tot, NGX, NGY, NGZ, grid_idx = parse_grid_info(content)
    pot = parse_potential(content, atom_tot, NGX, NGY, NGZ, grid_idx)

    if args.interp > 1:
        pot = interpolate_potential(pot, args.interp)
    logging.info(f"Final grid shape for plotting: {pot.shape}")

    plot_3d_potential(
        lattice, pot,
        sample=args.sample,
        vmin=args.vmin,
        vmax=args.vmax
    )

