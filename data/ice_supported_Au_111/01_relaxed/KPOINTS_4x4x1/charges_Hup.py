#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

z_cell = 49.2004013061999999

number_bilayer = [1,
                  2,2,
                  3,3,3,
                  4,4,4,4,
                  5,5,5,5,5,
                  6,6,6,6,6,6
                  ]

hydrogen_charges = [0.3380208958333333,
                    0.33343835416666673, 0.3408492291666667,
                    0.32631383333333336, 0.31998516666666665, 0.3307930208333334,
                    0.31018468749999994,0.31812287499999997,0.31526810416666673,0.3455277708333333,
                    0.3019710416666666, 0.31126841666666677, 0.3142170208333334, 0.3127371666666666, 0.33607179166666673,
                    0.30178687499999995,0.31139518750000006,0.31443189583333336,0.3130669791666667,0.3360273541666668,0.34,
                    ]

hydrogen_charges = [-x + 1 for x in hydrogen_charges]

oxygen_charges = [7.328148958333333,
                  7.3332990833333325, 7.355318708333333,
                  7.343306666666668, 7.359484999999999, 7.393085541666667,
                  7.369388583333333,7.3639823333333325,7.371821374999999,7.368475291666667,
                  7.386314541666667, 7.378184833333332, 7.371747166666665, 7.374628875, 7.3956021666666665,
                  7.3843305,7.377955208333333,7.371309041666668,7.3742492916666675,7.394824833333331,7.40,
                  ]

oxygen_charges = [-x + 6 for x in oxygen_charges]

hyd_plus_ox = [2*h + o for h,o in zip(hydrogen_charges,oxygen_charges)]

colors = [
    "#377eb8",
    "#377eb8","#e41a1c",
    "#377eb8","#e41a1c","#4daf4a",
    "#377eb8","#e41a1c","#4daf4a","#ff7f00",
    "#377eb8","#e41a1c","#4daf4a","#ff7f00","#984ea3",
    "#377eb8","#e41a1c","#4daf4a","#ff7f00","#984ea3","#a65628",
    ]

labels = ['Bilayer 1','Bilayer 2','Bilayer 3','Bilayer 4','Bilayer 5','Bilayer 6']

# Create the figure and subplots
fig, axs = plt.subplots(3, 1, figsize=(8,12))
fig.subplots_adjust(hspace=0.3)

# Plot Hydrogen Charges
axs[0].set_ylabel('H Charge (e)')
axs[0].xaxis.set_minor_locator(AutoMinorLocator())

for i in range(len(number_bilayer)):
    if i <= 14:
        axs[0].scatter(number_bilayer[i], hydrogen_charges[i], color=colors[i], s=100)
    else:
        axs[0].scatter(number_bilayer[i], hydrogen_charges[i], color=colors[i], s=100, label=labels[i - 15])

axs[0].legend(loc='upper right')

# Plot Oxygen Charges
axs[1].set_ylabel('O Charge (e)')
axs[1].xaxis.set_minor_locator(AutoMinorLocator())

for i in range(len(number_bilayer)):
    marker_style = 'v' if i % 2 == 0 else '^'  # Alternate between triangle markers
    axs[1].scatter(number_bilayer[i], oxygen_charges[i], color=colors[i], s=100)

# Plot H2O Charges (Hydrogen + Oxygen)
axs[2].set_xlabel('Number of Bilayers')
axs[2].set_ylabel('Charge (e)')

for i in range(len(number_bilayer)):
    marker_style = 'D' if i % 2 == 0 else '*'  # Alternate between diamond and star markers
    axs[2].scatter(number_bilayer[i], hyd_plus_ox[i], color=colors[i], s=100)

# Show plot
plt.show()
