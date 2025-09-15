#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

def read_com_potentials(filename="COM_potentials_summary.txt"):
    potentials = []
    with open(filename, "r") as f:
        for line in f:
            if "Potential at water COM" in line:
                try:
                    val_str = line.split(":")[-1].strip().split()[-2]
                    potentials.append(float(val_str))
                except Exception as e:
                    print(f"Skipping line (parse error): {line.strip()} — {e}")
    return np.array(potentials)

def plot_histogram(values):
    mean_val = np.mean(values)
    std_val = np.std(values)

    plt.figure(figsize=(7,5))

    # Histogram (counts)
    n_bins = 15
    counts, bins, patches = plt.hist(values, bins=n_bins, color='steelblue', edgecolor='black',
                                     alpha=0.6, label='Histogram (count)')

    # KDE rescaled to histogram counts
    kde = gaussian_kde(values)
    x_vals = np.linspace(min(bins), max(bins), 300)
    scale_factor = len(values) * (bins[1] - bins[0])
    plt.plot(x_vals, kde(x_vals) * scale_factor, color='darkred', linewidth=2, label="KDE")

    # Statistical indicators
    plt.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f"Mean = {mean_val:.2f} eV")
    plt.axvline(mean_val + std_val, color='orange', linestyle='--', linewidth=1.5, label=f"±σ = {std_val:.2f} eV")
    plt.axvline(mean_val - std_val, color='orange', linestyle='--', linewidth=1.5)

    # Labels
    plt.xlabel("Electrostatic potential at COM (eV)", fontsize=14)
    plt.ylabel("Count", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)

    stats_text = f"Mean = {mean_val:.3f} eV\nStd Dev = {std_val:.3f} eV"
    plt.gca().text(0.95, 0.95, stats_text, ha='right', va='top',
                   transform=plt.gca().transAxes,
                   fontsize=12,
                   bbox=dict(boxstyle="round", facecolor='white', edgecolor='gray'))

    plt.savefig("histogram_com_potential.png", dpi=300)
    plt.show()

def main():
    values = read_com_potentials()
    if len(values) == 0:
        print("No potential values found.")
        return
    plot_histogram(values)

if __name__ == "__main__":
    main()

