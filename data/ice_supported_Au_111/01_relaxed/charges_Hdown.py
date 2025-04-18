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

hydrogen_charges = [0.3278499999999999,
                    0.3354755208333333,0.31776197916666665,
                    0.3315077083333334,0.3176482708333333,0.3077466458333333,
                    0.323266625,0.3182875833333333,0.31865237500000004,0.31729674999999996,
                    0.32022764583333335,0.32062900000000005,0.3118600416666667,0.3202412083333334,0.3101538333333333,
                    0.32316945833333327,0.3230257916666667,0.31548233333333325,0.3140941874999999,0.31678493750000003,0.3119736250000001,
                    ]

hydrogen_charges = [-x + 1 for x in hydrogen_charges]

oxygen_charges = [7.344255124999999,
                  7.333414000000001,7.3311261666666665,
                  7.344693374999999,7.363773958333333,7.3324074999999995,
                  7.360341416666668,7.364061666666668,7.362225375000001,7.304898208333333,
                  7.368341000000001,7.359670416666666,7.375425375,7.359960458333334,7.311622458333333,
                  7.361862583333335,7.355631125000002,7.369006041666668,7.3718475833333335,7.365984208333333,7.3054760833333345,
                  ]

oxygen_charges = [-x + 6 for x in oxygen_charges]

hyd_plus_ox = [2*h + o for h,o in zip(hydrogen_charges,oxygen_charges)]

colors = ['blue',
          'blue','red',
          'blue','red','green',
          'blue','red','green','orange',
          'blue','red','green','orange','purple',
          'blue','red','green','orange','purple','pink',
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
axs[2].set_ylabel('H$_2$O Charge (e)')
axs[2].xaxis.set_minor_locator(AutoMinorLocator())

for i in range(len(number_bilayer)):
    marker_style = 'D' if i % 2 == 0 else '*'  # Alternate between diamond and star markers
    axs[2].scatter(number_bilayer[i], hyd_plus_ox[i], color=colors[i], s=100)

# Show plot
plt.show()
