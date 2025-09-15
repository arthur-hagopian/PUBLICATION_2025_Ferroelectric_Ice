#!/usr/bin/env python3
import os, sys, logging
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# === USER SETTINGS: just edit these lines! ===
LOCPOT_FILE   = "LOCPOT"        # your LOCPOT file
INTERP_FACTOR = 1               # up-sampling (1 = no interp)
P1_FRAC = (0.6018553423875943, 0.9105603724948779, 0.0982981996760017)  # O
P2_FRAC = (0.5574453127387926, 0.0099732040061511, 0.0894567443376886)  # H1
P3_FRAC = (0.6897665793844041, 0.9184035497050402, 0.1251760079063473)  # H2
PLANE_MARGIN  = 2               # extend triangle by this fraction
GRID_RES      = 200             # resolution of the plane grid
VMIN, VMAX    = -5, 5           # color clamp (None = auto)
COORD_MODE    = "fractional"    # "fractional" or "cartesian"
# ===============================================

def read_locpot(fn):
    if not os.path.isfile(fn):
        logging.error(f"File {fn!r} not found."); sys.exit(1)
    lines = open(fn).readlines()
    logging.info("LOCPOT read.")
    return lines

def parse_lattice(lines):
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[2+i].split()] for i in range(3)]) * scale
    logging.info(f"Lattice parsed (scale={scale}).")
    return lat

def parse_grid(lines):
    counts = [int(x) for x in lines[6].split()]
    atom_tot = sum(counts)
    idx = 9 + atom_tot
    NGX, NGY, NGZ = map(int, lines[idx].split())
    logging.info(f"Grid dims: {NGX}×{NGY}×{NGZ}")
    return NGX, NGY, NGZ, idx

def parse_pot(lines, NGX, NGY, NGZ, idx):
    total = NGX*NGY*NGZ
    vals = []
    i = idx
    while len(vals) < total:
        i += 1
        vals += [float(x) for x in lines[i].split()]
    pot = np.zeros((NGX,NGY,NGZ), float)
    l = 0
    for k in range(NGZ):
        for j in range(NGY):
            for i in range(NGX):
                pot[i,j,k] = vals[l]; l+=1
    logging.info("Potential grid parsed.")
    return pot

def interpolate(pot, factor):
    if factor<=1: return pot
    NGX,NGY,NGZ = pot.shape
    x = np.arange(NGX)/NGX; y = np.arange(NGY)/NGY; z = np.arange(NGZ)/NGZ
    interp = RegularGridInterpolator((x,y,z), pot, method='linear',
                                     bounds_error=False, fill_value=None)
    new = [np.linspace(0,1,NGX*factor,endpoint=False),
           np.linspace(0,1,NGY*factor,endpoint=False),
           np.linspace(0,1,NGZ*factor,endpoint=False)]
    ii,jj,kk = np.meshgrid(*new, indexing='ij')
    pts = np.vstack([ii.ravel(), jj.ravel(), kk.ravel()]).T
    vals = interp(pts).reshape(NGX*factor,NGY*factor,NGZ*factor)
    logging.info(f"Up-sampled to {vals.shape}.")
    return vals

def build_interp_func(pot):
    NGX,NGY,NGZ = pot.shape
    x = np.arange(NGX)/NGX; y = np.arange(NGY)/NGY; z = np.arange(NGZ)/NGZ
    return RegularGridInterpolator((x,y,z), pot,
                                  method='linear',
                                  bounds_error=False,
                                  fill_value=None)

def plot_plane(lat, pot):
    def circle_mask(center, radius, xgrid, ygrid):
        return (xgrid - center[0])**2 + (ygrid - center[1])**2 <= radius**2

    p1 = np.array(P1_FRAC)
    p2 = np.array(P2_FRAC)
    p3 = np.array(P3_FRAC)

    # build in-plane basis in fractional space
    d1 = p2 - p1
    d2 = p3 - p1
    M  = np.column_stack((d1,d2))
    uv = np.linalg.lstsq(M, (np.vstack((p1,p2,p3)).T - p1[:,None]), rcond=None)[0]
    u_coords, v_coords = uv[0], uv[1]

    # find extents and add margin
    umin, umax = u_coords.min(), u_coords.max()
    vmin_, vmax_ = v_coords.min(), v_coords.max()
    du = (umax-umin)*PLANE_MARGIN
    dv = (vmax_-vmin_)*PLANE_MARGIN
    u_range = np.linspace(umin-du, umax+du, GRID_RES)
    v_range = np.linspace(vmin_-dv, vmax_+dv, GRID_RES)
    uu, vv  = np.meshgrid(u_range, v_range, indexing='ij')

    # build every frac-space point on plane: p(u,v) = p1 + u*d1 + v*d2
    pts_frac = p1 + uu[...,None]*d1 + vv[...,None]*d2
    pts_flat = pts_frac.reshape(-1,3)

    # build and use the interpolator
    f = build_interp_func(pot)

    # Compute center of mass (fractional coordinates)
    m_O, m_H = 15.999, 1.008
    total_mass = m_O + 2 * m_H
    com_frac = ((p1 * m_O + p2 * m_H + p3 * m_H) / total_mass) % 1.0

    pot_COM = float(f(com_frac)[0])
    logging.info(f"Potential at water COM (O+H1+H2) : {pot_COM:.6f} eV")

    # compute potential map
    vals = f(pts_flat).reshape(GRID_RES, GRID_RES)

    # axes choice
    if COORD_MODE == "cartesian":
        cart = pts_flat @ lat
        xs = cart[:,0].reshape(GRID_RES, GRID_RES)
        ys = cart[:,1].reshape(GRID_RES, GRID_RES)
        extent = [xs.min(), xs.max(), ys.min(), ys.max()]
        xlabel, ylabel = "X (Å)", "Y (Å)"
        coords = [(p1 @ lat)[:2], (p2 @ lat)[:2], (p3 @ lat)[:2], (com_frac @ lat)[:2]]
    else:
        extent = [u_range.min(), u_range.max(), v_range.min(), v_range.max()]
        xlabel, ylabel = "u (fractional)", "v (fractional)"
        coords = [
            (u_coords[0], v_coords[0]),
            (u_coords[1], v_coords[1]),
            (u_coords[2], v_coords[2]),
            (np.nan, np.nan)  # COM will be computed below
        ]
        # Solve COM projection
        com_vec = com_frac - p1
        coeffs = np.linalg.lstsq(np.column_stack((d1,d2)), com_vec, rcond=None)[0]
        coords[3] = (coeffs[0], coeffs[1])

    labels = ["O", "H1", "H2", "COM"]

    # plot heatmap
    fig, ax = plt.subplots(figsize=(6,5))
    im = ax.imshow(
        vals.T, origin='lower', extent=extent,
        cmap='seismic', vmin=VMIN, vmax=VMAX, aspect='equal'
    )
    fig.colorbar(im, ax=ax, label='Potential (eV)')
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title("LOCPOT on plane through O, H1, H2")

    for (x,y,label) in zip(coords, coords, labels):
        ax.plot(x[0], x[1], marker='+', color='k', markersize=12, mew=2)
        ax.text(x[0], x[1], " "+label, color='k', fontsize=12,
                verticalalignment='bottom', horizontalalignment='left')

    plt.tight_layout()
    plt.show()

def main():
    lines = read_locpot(LOCPOT_FILE)
    lat   = parse_lattice(lines)
    NGX,NGY,NGZ,idx = parse_grid(lines)
    pot   = parse_pot(lines, NGX,NGY,NGZ, idx)
    if INTERP_FACTOR>1:
        pot = interpolate(pot, INTERP_FACTOR)
    logging.info(f"Grid shape now {pot.shape}")
    plot_plane(lat, pot)

if __name__=='__main__':
    main()

