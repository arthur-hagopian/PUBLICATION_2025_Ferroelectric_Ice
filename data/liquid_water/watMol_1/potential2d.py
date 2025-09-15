#!/usr/bin/env python3
import os, sys, logging
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# === USER SETTINGS: just edit these lines! ===
LOCPOT_FILE   = "LOCPOT"        # your LOCPOT file
INTERP_FACTOR = 1               # up-sampling (1 = no interp)
P1_FRAC       = (0.5800782483774682, 0.3846091147133421, 0.1304954711713791)  # O
P2_FRAC       = (0.5339218942628037, 0.4701148254184986, 0.1217281887742132)  # H1
P3_FRAC       = (0.6784725264648657, 0.4068722045053974, 0.0818958929790412)  # H2
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
    from matplotlib.path import Path

    def circle_mask(center, radius, xgrid, ygrid):
        return (xgrid - center[0])**2 + (ygrid - center[1])**2 <= radius**2

    # fractional “atom” points
    p1 = np.array(P1_FRAC)
    p2 = np.array(P2_FRAC)
    p3 = np.array(P3_FRAC)

    # build in-plane basis in fractional space
    d1 = p2 - p1
    d2 = p3 - p1
    M  = np.column_stack((d1,d2))  # 3×2
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

    # ——— LOG POTENTIALS AT P1, P2, P3 AND MIDPOINT ———
    pot_O  = float(f(p1))
    pot_H1 = float(f(p2))
    pot_H2 = float(f(p3))
    mid_frac = ((p2 + p3) / 2.0) % 1.0
    pot_mid = float(f(mid_frac))
    diff    = pot_mid - pot_O

    logging.info(f"Potential at O  (P1)           : {pot_O:.6f} eV")
    logging.info(f"Potential at H1 (P2)           : {pot_H1:.6f} eV")
    logging.info(f"Potential at H2 (P3)           : {pot_H2:.6f} eV")
    logging.info(f"Potential at midpoint (H1–H2)  : {pot_mid:.6f} eV")
    logging.info(f"Difference (midpoint – O)      : {diff:.6f} eV")
    # ———————————————————————————————

    # compute potential map
    vals = f(pts_flat).reshape(GRID_RES, GRID_RES)

    # axes choice
    if COORD_MODE == "cartesian":
        cart = pts_flat @ lat
        xs = cart[:,0].reshape(GRID_RES, GRID_RES)
        ys = cart[:,1].reshape(GRID_RES, GRID_RES)
        extent = [xs.min(), xs.max(), ys.min(), ys.max()]
        xlabel, ylabel = "X (Å)", "Y (Å)"
        p1_xy = (p1 @ lat)[:2]
        p2_xy = (p2 @ lat)[:2]
        p3_xy = (p3 @ lat)[:2]
        circle_radius = 1.75  # in Angstrom
        xgrid, ygrid = xs, ys
    else:
        extent = [u_range.min(), u_range.max(), v_range.min(), v_range.max()]
        xlabel, ylabel = "u (fractional)", "v (fractional)"
        p1_xy = (u_coords[0], v_coords[0])
        p2_xy = (u_coords[1], v_coords[1])
        p3_xy = (u_coords[2], v_coords[2])
        circle_radius = 1.75 / np.linalg.norm(d1 @ lat[:,:2])  # rough approx
        xgrid, ygrid = uu, vv

    # generate union of circles mask
    mask1 = circle_mask(p1_xy, circle_radius, xgrid, ygrid)
    mask2 = circle_mask(p2_xy, circle_radius, xgrid, ygrid)
    mask3 = circle_mask(p3_xy, circle_radius, xgrid, ygrid)
    union_mask = mask1 | mask2 | mask3

    # plot heatmap
    fig, ax = plt.subplots(figsize=(6,5))
    im = ax.imshow(
        vals.T, origin='lower', extent=extent,
        cmap='seismic', vmin=VMIN, vmax=VMAX, aspect='equal'
    )
    fig.colorbar(im, ax=ax, label='Potential (eV)')
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title("LOCPOT on plane through P1,P2,P3")

    # draw dotted contour of the union
    union_mask_T = union_mask.T.astype(int)
    ax.contour(xgrid, ygrid, union_mask_T, levels=[0.5],
               colors='k', linewidths=2, linestyles='dotted')

    # annotate atoms
    for (x,y,label) in [
        (p1_xy[0], p1_xy[1], "O"),
        (p2_xy[0], p2_xy[1], "H1"),
        (p3_xy[0], p3_xy[1], "H2")
    ]:
        ax.plot(x, y, marker='+', color='k', markersize=12, mew=2)
        ax.text(x, y, " "+label, color='k', fontsize=12,
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

