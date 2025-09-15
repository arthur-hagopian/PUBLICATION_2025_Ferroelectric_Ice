#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# Constant 
z_cell = 49.2004013061999999

# Charges
number_bilayer = [1,
                  2,2,
                  3,3,3,
                  4,4,4,4,
                  5,5,5,5,5,
                  6,6,6,6,6,6
                  ]

hydrogen_charges = [0.3065033333333333,
                    0.30795389583333327,0.33181366666666673,
                    0.30525193749999996,0.3108988124999999,0.3107817708333333,
                    0.30725895833333344,0.30949941666666675,0.3086407291666667,0.3094771666666667,
                    0.3086270833333334,0.4988082083333336,0.33069033333333336,0.3109645625,0.31108812499999994,
                    0.30798937500000007,0.30733035416666665,0.2792874791666667,0.3057228333333334,0.3059234374999999,0.31743789583333326,
                    ]

hydrogen_charges = [-x + 1 for x in hydrogen_charges]

oxygen_charges = [7.386992416666668,
                  7.397338375000001,7.323125750000002,
                  7.410381333333334,7.378305166666668,7.357452458333332,
                  7.4169633333333325,7.3814615833333335,7.379793208333335,7.352034291666668,
                  7.458164750000002,6.911291208333334,7.372961916666668,7.379155208333334,7.358075666666667,
                  7.424469124999999,7.3864394583333315,7.441814916666668,7.389002958333332,7.3839939166666655,7.326905374999999,
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
fig, axs = plt.subplots(3, 1)
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
axs[2].set_ylim(-0.15,0.15)
axs[2].set_xlabel('Number of Bilayers')
axs[2].set_ylabel('Charge (e)')
#axs[2].xaxis.set_minor_locator(AutoMinorLocator())

for i in range(len(number_bilayer)):
    marker_style = 'D' if i % 2 == 0 else '*'  # Alternate between diamond and star markers
    axs[2].scatter(number_bilayer[i], hyd_plus_ox[i], color=colors[i], s=100)

# Show plot
plt.show()
