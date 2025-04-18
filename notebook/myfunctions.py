# MODULES
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import plotly.graph_objects as go
from pathlib import Path
import logging
import sys

### LOGGING

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def configure_logger(level=logging.INFO):
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)

    logger = logging.getLogger()
    logger.handlers = [handler]  # Replace existing handlers
    logger.setLevel(level)

### READ & LOAD DATA

def read_contcar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    header = lines[:9]
    
    # Extract element symbols and number of atoms per element
    element_symbols = header[5].split()
    element_counts = list(map(int, header[6].split()))
    
    total_atoms = sum(element_counts)
    start_index = 9
    end_index = start_index + total_atoms

    # Read atomic positions
    atomic_positions = np.array([
        list(map(float, line.split()[:3])) 
        for line in lines[start_index:end_index]
    ])

    # Map each element to its slice of positions
    element_positions = {}
    idx = 0
    for symbol, count in zip(element_symbols, element_counts):
        element_positions[symbol] = atomic_positions[idx:idx+count]
        idx += count

    return element_positions

def load_positions(file_paths):
    H_result = {}
    O_result = {}
    Au_result = {}
    for path in file_paths:
        positions = read_contcar(path)
        H_positions = positions.get("H")
        O_positions = positions.get("O")
        Au_positions = positions.get("Au")
        H_result[path] = H_positions
        O_result[path] = O_positions
        Au_result[path] = Au_positions

    return H_result, O_result, Au_result

def read_doscar(path):
    with open(path, 'r') as f:
        lines = f.readlines()
        header = lines[:6]
        num_atoms = int(header[0].split()[0])
        efermi = float(header[5].split()[3])
        ndos = int(header[5].split()[2])
        total_dos_data = np.array([list(map(float,line.split())) for line in lines[6:6+ndos]])
        energy = total_dos_data[:,0] - efermi
        total_dos = total_dos_data[:,1]
        pdos = []
        start_idx = 6 + ndos + 1
        for atom in range(num_atoms):
            atom_dos_data = []
            for line in lines[start_idx:start_idx+ndos]:
                columns = list(map(float, line.split()))
                atom_dos_data.append(columns[1:])
            pdos.append(np.array(atom_dos_data))
            start_idx += ndos + 1

    return {
        "energy": energy,
        "pdos": pdos
    }

def load_doscar(file_paths):
    result = {}
    for path in file_paths:
        value = read_doscar(path)
        if value is not None:
            result[path] = value
    return result

def read_bader(filename):
    charges = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # skip blank lines, comment‐headers, or separator lines
            if not line or line.startswith('#') or all(c in '- ' for c in line):
                continue

            cols = line.split()
            # we expect at least 5 columns: index, X, Y, Z, CHARGE, ...
            if len(cols) < 5:
                logger.debug(f"Skipping malformed line: {line}")
                continue

            # first column should be atom index
            try:
                _ = int(cols[0])
            except ValueError:
                logger.debug(f"Skipping non‐data line: {line}")
                continue

            # parse the CHARGE (5th column, zero‐based index 4)
            try:
                charge = float(cols[4])
            except ValueError:
                logger.warning(f"Could not parse charge on line: {line}")
                continue

            charges.append(charge)

    return charges

def load_bader(file_paths):
    result = {}
    for path in file_paths:
        value = read_bader(path)
        if value is not None:
            result[path] = value
    return result

def collect_data_files(base_path):
    data_files = {
        "CONTCAR": [],
        "DOSCAR": [],
        "ACF.dat": [],
    }
    for file_path in Path(base_path).rglob("*"):
        if file_path.name in data_files:
            data_files[file_path.name].append(file_path)
    return data_files

def load_all_data(base_path):
    
    files = collect_data_files(base_path)

    H_positions, O_positions, Au_positions = load_positions(files["CONTCAR"])
    pdos = load_doscar(files["DOSCAR"])
    bader = load_bader(files["ACF.dat"])

    return {
        "H_positions": H_positions,
        "O_positions": O_positions,
        "Au_positions": Au_positions,
        "pdos": pdos,
        "bader": bader,
    }

### COMPUTE AND PLOT

def slice_z_layers(h_positions, o_positions, au_positions, z_bounds):

    base_colors = ['blue','red','green','orange','purple','pink']
    base_alphas = [1.0, 0.5]  # H, O

    n_au = len(au_positions) if au_positions is not None else 0
    n_h = len(h_positions)
    groups, labels, colors, alphas = [], [], [], []

    n_slabs = len(z_bounds) - 1
    for i in range(n_slabs):
        z_lo, z_hi = z_bounds[i], z_bounds[i+1]

        # hydrogen indices in this slab
        h_idx = [
            n_au+idx+1
            for idx, pos in enumerate(h_positions)
            if z_lo < pos[2] <= z_hi
        ]
        # oxygen indices (offset by n_h)
        o_idx = [
            n_au + n_h + idx+1
            for idx, pos in enumerate(o_positions)
            if z_lo < pos[2] <= z_hi
        ]

        groups += [h_idx, o_idx]
        labels += [f"Bilayer {i+1} : H", f"Bilayer {i+1} : O"]
        colors += [base_colors[i % len(base_colors)]] * 2
        alphas += base_alphas

    return groups, labels, colors, alphas

def average_pdos_by_group(pdos, groups):
    """
    For each group (list of atom‐indices), sum & average the
    PDOS rows over atoms, returning a list of 1D arrays.
    """
    averaged = []
    for group in groups:
        if not group:
            averaged.append(None)
            continue
        # sum over atoms, remembering pdos is zero‑based
        stack = [pdos[idx-1].sum(axis=1) for idx in group]
        avg = np.mean(stack, axis=0)
        averaged.append(avg)
    return averaged

def plot_pdos(energy, averaged_dos, labels, colors, alphas,
             figsize=(5,3), xlim=(-15,5), ylim_top=3.5):
    fig, ax = plt.subplots(figsize=figsize)
    for y, lbl, c, a in zip(averaged_dos, labels, colors, alphas):
        if y is None:
            continue
        ax.plot(energy, y, label=lbl, color=c, alpha=a)

    ax.set_xlabel(r"E - E$_{F}$ (eV)")
    ax.set_ylabel("PDOS (a.u.)")
    ax.set_xlim(xlim)
    ax.set_ylim(top=ylim_top)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax.legend()
    plt.tight_layout()
    plt.show()

def plot_pdos_all_bilayers(data, base_path, folder_label, z_bounds, max_layers=6):
    for n in range(1, max_layers+1):
        logger.info(f"--- bilayer_{n} ---")
        # find the matching DOSCAR key & load positions
        dos_key = next(
            k for k in data["pdos"]
            if f"{base_path}/{folder_label}_{n}" in str(k)
        )
        dos = data["pdos"][dos_key]
        h_pos = data["H_positions"].get(dos_key.parent/"CONTCAR")
        o_pos = data["O_positions"].get(dos_key.parent/"CONTCAR")
        au_pos = data["Au_positions"].get(dos_key.parent/"CONTCAR")

        if h_pos is None or o_pos is None:
            logger.warning("Missing positions – skipping")
            continue

        # 1) select
        groups, labels, colors, alphas = slice_z_layers(h_pos, o_pos, au_pos, z_bounds)
        # 2) process
        avg_dos = average_pdos_by_group(dos["pdos"], groups)
        # 3) plot
        plot_pdos(dos["energy"], avg_dos, labels, colors, alphas)

def plot_bilayer_avg_charge(data, base_path, folder_label, z_bounds, max_layers=6):
    """

    """
    base_colors = ['blue', 'red', 'green', 'orange', 'purple', 'pink']
    fig, ax = plt.subplots(figsize=(6, 4))

    # Track which layer labels have been added to the legend
    legend_done = set()

    for n in range(1, max_layers + 1):
        # locate Bader file for this bilayer_n
        bader_key = next(
            (k for k in data['bader'] if f"{base_path}/{folder_label}_{n}" in str(k.parent)),
            None
        )
        if bader_key is None:
            logger.warning(f"No Bader data for bilayer_{n}, skipping")
            continue

        charges = data['bader'][bader_key]
        contcar = bader_key.parent / 'CONTCAR'
        h_pos = data['H_positions'].get(contcar)
        o_pos = data['O_positions'].get(contcar)
        au_pos = data['Au_positions'].get(contcar)
        if h_pos is None or o_pos is None:
            logger.warning(f"Missing positions for bilayer_{n}, skipping")
            continue

        # split into n layers
        groups, labels, colors, alphas = slice_z_layers(h_pos, o_pos, au_pos, z_bounds)
        # groups are [H_layer1, O_layer1, H_layer2, O_layer2, ...]

        for layer in range(n):
            # combine H and O indices for this layer
            idx_h = groups[2*layer]
            idx_o = groups[2*layer + 1]
            avg_charge_h = np.mean([- charges[i-1] + 1 for i in idx_h])
            avg_charge_o = np.mean([- charges[i-1] + 6 for i in idx_o])

            avg_charge = 2*avg_charge_h + avg_charge_o

            # scatter at x = n
            c = colors[2*layer]  # same color for both H/O of layer
            label = f"Bilayer {layer+1}"
            if label not in legend_done:
                ax.scatter(n, avg_charge, color=c, label=label, zorder=3)
                legend_done.add(label)
            else:
                ax.scatter(n, avg_charge, color=c, zorder=3)

    ax.set_xlabel('Number of bilayers (n)')
    ax.set_ylabel('Average Bader charge per bilayer (e)')   
    ax.set_xticks(range(1, max_layers + 1))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.legend()
    ax.grid()
    plt.tight_layout()
    plt.show()

def plot_interbilayer_distance(data, base_path, z_bounds, lz, max_layers=6):
    """

    """
    base_colors = ['blue', 'red', 'green', 'orange', 'purple', 'pink']
    fig, ax = plt.subplots(figsize=(6, 4))
    legend_done = set()

    for n in range(2, max_layers + 1):
        bader_key = next(
            (k for k in data['bader'] if f"{base_path}/bilayer_{n}" in str(k.parent)),
            None
        )
        if bader_key is None:
            logger.warning(f"No Bader data for bilayer_{n}, skipping")
            continue

        contcar = bader_key.parent / 'CONTCAR'
        h_pos = data['H_positions'].get(contcar)
        o_pos = data['O_positions'].get(contcar)
        au_pos = data['Au_positions'].get(contcar)

        if h_pos is None or o_pos is None:
            logger.warning(f"Missing positions for bilayer_{n}, skipping")
            continue

        # Ensure arrays have correct shape
        au_pos = np.array(au_pos).reshape(-1, 3) if au_pos is not None else np.empty((0, 3))
        h_pos = np.array(h_pos).reshape(-1, 3)
        o_pos = np.array(o_pos).reshape(-1, 3)

        # Layer assignment
        groups, labels, colors, alphas = slice_z_layers(h_pos, o_pos, au_pos, z_bounds)
        pos_all = np.concatenate([au_pos, h_pos, o_pos])

        # Compute z positions (converted to angstrom)
        layer_z_avgs = []
        for layer in range(n):
            idx_h = [i - 1 for i in groups[2 * layer]]
            idx_o = [i - 1 for i in groups[2 * layer + 1]]
            if not idx_h or not idx_o:
                logger.warning(f"Empty group in bilayer_{n}, layer {layer+1}, skipping")
                continue
            try:
                z_h = pos_all[idx_h, 2] * lz
                z_o = pos_all[idx_o, 2] * lz
                avg_z = np.mean(np.concatenate([z_h, z_o]))
                layer_z_avgs.append(avg_z)
            except IndexError as e:
                logger.error(f"Index error for bilayer_{n}, layer {layer+1}: {e}")
                continue

        if len(layer_z_avgs) < 2:
            logger.warning(f"Not enough valid layers in bilayer_{n}, skipping")
            continue

        layer_z_avgs = sorted(layer_z_avgs)
        interlayer_distances = np.diff(layer_z_avgs)

        for i, d in enumerate(interlayer_distances):
            c = base_colors[i % len(base_colors)]
            label = f"Bilayer {i+1}"
            if label not in legend_done:
                ax.scatter(n, d, color=c, label=label, zorder=3)
                legend_done.add(label)
            else:
                ax.scatter(n, d, color=c, zorder=3)

    ax.set_xlabel('Number of bilayers (n)')
    ax.set_ylabel('Interbilayer distance (Å)')
    ax.set_xticks(range(2, max_layers + 1))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.legend()
    ax.grid()
    plt.tight_layout()
    plt.show()






# VISUALIZE

def plotly_3d_atoms(H_pos, O_pos, Au_pos=None, cell=None, projection='perspective', eye=(1.5, 1.5, 1.0)):
    data = []

    if H_pos is not None and len(H_pos) > 0:
        data.append(go.Scatter3d(
            x=H_pos[:, 0], y=H_pos[:, 1], z=H_pos[:, 2],
            mode='markers',
            marker=dict(
                size=4, color='white',
                line=dict(color='black', width=1)
            ),
            name='Hydrogen'
        ))

    if O_pos is not None and len(O_pos) > 0:
        data.append(go.Scatter3d(
            x=O_pos[:, 0], y=O_pos[:, 1], z=O_pos[:, 2],
            mode='markers',
            marker=dict(
                size=5, color='red',
                line=dict(color='black', width=1)
            ),
            name='Oxygen'
        ))

    if Au_pos is not None and len(Au_pos) > 0:
        data.append(go.Scatter3d(
            x=Au_pos[:, 0], y=Au_pos[:, 1], z=Au_pos[:, 2],
            mode='markers',
            marker=dict(
                size=6, color='gold',
                line=dict(color='black', width=1)
            ),
            name='Gold'
        ))

    # Cell box drawing
    if cell is not None:
        corners = np.array([
            [0, 0, 0], cell[0], cell[1], cell[2],
            cell[0] + cell[1], cell[0] + cell[2], cell[1] + cell[2],
            cell[0] + cell[1] + cell[2],
        ])
        edge_indices = [
            (0, 1), (0, 2), (0, 3), (1, 4), (1, 5),
            (2, 4), (2, 6), (3, 5), (3, 6),
            (4, 7), (5, 7), (6, 7)
        ]
        for i, j in edge_indices:
            data.append(go.Scatter3d(
                x=[corners[i][0], corners[j][0]],
                y=[corners[i][1], corners[j][1]],
                z=[corners[i][2], corners[j][2]],
                mode='lines',
                line=dict(color='gray', width=2),
                showlegend=False
            ))
        min_xyz = np.min(corners, axis=0)
        max_xyz = np.max(corners, axis=0)
    else:
        all_points = [arr for arr in [H_pos, O_pos, Au_pos] if arr is not None and len(arr) > 0]
        all_points = np.concatenate(all_points)
        min_xyz = np.min(all_points, axis=0)
        max_xyz = np.max(all_points, axis=0)

    fig = go.Figure(data=data)
    fig.update_layout(
        scene=dict(
            xaxis=dict(title='X (Å)', range=[min_xyz[0], max_xyz[0]]),
            yaxis=dict(title='Y (Å)', range=[min_xyz[1], max_xyz[1]]),
            zaxis=dict(title='Z (Å)', range=[min_xyz[2], max_xyz[2]]),
            aspectmode='data',
            camera=dict(
                projection=dict(type=projection),
                eye=dict(x=eye[0], y=eye[1], z=eye[2])
            )
        ),
        margin=dict(l=0, r=0, b=0, t=30)
    )
    fig.show()


