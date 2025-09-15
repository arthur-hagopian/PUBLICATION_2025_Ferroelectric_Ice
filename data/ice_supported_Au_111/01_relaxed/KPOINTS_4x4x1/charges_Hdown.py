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

hydrogen_charges = [0.3281387916666667,
                    0.33571235416666667, 0.3168580625,
                    0.3318309375, 0.3173590208333332, 0.3074921041666666,
                    0.3233780833333334, 0.31815131250000006, 0.3180254166666666, 0.3169662499999999,
                    0.3206307916666666, 0.32066045833333334, 0.31205395833333327, 0.3204937291666667, 0.3108041666666666,
                    0.3214695625, 0.32289518750000007, 0.3142047916666667, 0.3140418333333333, 0.31665672916666665, 0.31191597916666663,
                    ]

hydrogen_charges = [-x + 1 for x in hydrogen_charges]

oxygen_charges = [7.345594874999997,
                  7.334994208333334, 7.3304020833333325,
                  7.345882583333332, 7.364343708333333, 7.3313935,
                  7.361929166666666, 7.36433925, 7.363459125000002, 7.304452083333332,
                  7.367576458333335, 7.359625083333334, 7.375002583333335, 7.359784208333334, 7.311344541666668,
                  7.366922291666666, 7.355920625, 7.371522625000001, 7.371977208333335, 7.3662133333333335, 7.304889958333335,
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
