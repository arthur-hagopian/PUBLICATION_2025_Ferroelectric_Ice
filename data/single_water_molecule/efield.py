import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')


def read_locpot(filename):
    if not os.path.isfile(filename):
        logging.error(f"File {filename} not found.")
        sys.exit(1)

    with open(filename) as f:
        content = f.readlines()
    logging.info("LOCPOT file read successfully.")
    return content


def parse_lattice(content):
    scaling_factor = float(content[1].split()[0])
    lattice = np.zeros((3, 3))
    for i in range(3):
        lattice[i] = [float(x) for x in content[2 + i].split()]
    lattice *= scaling_factor

    norms = np.linalg.norm(lattice, axis=1)
    logging.info("Lattice vectors and scaling factor parsed.")
    return lattice, norms


def parse_grid_info(content):
    atom_counts = [int(x) for x in content[6].split()]
    atom_tot = sum(atom_counts)

    NGX, NGY, NGZ = [int(x) for x in content[9 + atom_tot].split()]
    logging.info(f"Grid size parsed: NGX={NGX}, NGY={NGY}, NGZ={NGZ}")
    return atom_tot, NGX, NGY, NGZ


def parse_potential(content, atom_tot, NGX, NGY, NGZ):
    potential = []
    i = 9 + atom_tot
    grid_point_tot = NGX * NGY * NGZ

    while len(potential) < grid_point_tot:
        i += 1
        potential += [float(x) for x in content[i].split()]

    l = 0
    potential_grid = np.zeros(shape=(NGX,NGY,NGZ))
    for k in range(NGZ):
        for j in range(NGY):
            for i in range(NGX):
                potential_grid[i,j,k] = potential[l]
                l += 1

    logging.info("Potential data parsed and reshaped.")
    return potential_grid


def compute_electric_field(potential_grid, resolutions):
    grad_x, grad_y, grad_z = np.gradient(potential_grid, *resolutions)
    grad_magnitude = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
    logging.info("Electric field magnitude computed.")
    return grad_magnitude


def extract_plane(grid, plane_direction, plane_fraction):
    NGX, NGY, NGZ = grid.shape
    
    if plane_direction == 'x':
        ix = int(plane_fraction * NGX)
        plane = grid[ix, :, :]
    elif plane_direction == 'y':
        iy = int(plane_fraction * NGY)
        plane = grid[:, iy, :]
    elif plane_direction == 'z':
        iz = int(plane_fraction * NGZ)
        plane = grid[:, :, iz]
    else:
        logging.error("Invalid plane direction. Choose 'x', 'y', or 'z'.")
        sys.exit(1)

    logging.info(f"Plane extracted along {plane_direction}-axis at fraction {plane_fraction}.")
    return plane


def plot_planes(potential_plane, grad_plane, axis_lengths, title="Potential and Electric Field Slices"):
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))

    extent = [0, axis_lengths[0], 0, axis_lengths[1]]

    cs1 = axs[0].imshow(potential_plane.T, origin='lower', extent=extent, aspect='auto', cmap="RdYlBu_r", interpolation='bicubic')
    axs[0].set_title("Electrostatic Potential")
    axs[0].set_xlabel("X (Å)")
    axs[0].set_ylabel("Y (Å)")
    axs[0].set_xlim(13,17)
    axs[0].set_ylim(13,16.5)
    fig.colorbar(cs1, ax=axs[0])

    cs2 = axs[1].imshow(grad_plane.T, origin='lower', extent=extent, aspect='auto', cmap="RdYlBu_r", interpolation='bicubic')
    axs[1].set_title("Electric Field")
    axs[1].set_xlabel("X (Å)")
    axs[1].set_ylabel("Y (Å)")
    axs[1].set_xlim(13,17)
    axs[1].set_ylim(13,16.5)
    fig.colorbar(cs2, ax=axs[1])

    plt.tight_layout()
    plt.show()


def main(filename, plane_direction, plane_fraction):
    content = read_locpot(filename)
    lattice, norms = parse_lattice(content)
    atom_tot, NGX, NGY, NGZ = parse_grid_info(content)
    resolutions = norms[0]/NGX, norms[1]/NGY, norms[2]/NGZ
    potential_grid = parse_potential(content, atom_tot, NGX, NGY, NGZ)
    grad_magnitude = compute_electric_field(potential_grid, resolutions)

    potential_plane = extract_plane(potential_grid, plane_direction, plane_fraction)
    grad_plane = extract_plane(grad_magnitude, plane_direction, plane_fraction)

    if plane_direction == 'x':
        axis_lengths = (norms[1], norms[2])
    elif plane_direction == 'y':
        axis_lengths = (norms[0], norms[2])
    elif plane_direction == 'z':
        axis_lengths = (norms[0], norms[1])

    plot_planes(potential_plane, grad_plane, axis_lengths, title=f"Slices at {plane_direction}={plane_fraction}")


if __name__ == "__main__":
    plane_direction = 'z'  # choose between 'x', 'y', 'z'
    plane_fraction = 0.5   # between 0 and 1

    main(filename="LOCPOT", plane_direction=plane_direction, plane_fraction=plane_fraction)
