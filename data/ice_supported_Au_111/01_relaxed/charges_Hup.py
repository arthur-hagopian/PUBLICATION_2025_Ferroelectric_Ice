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

hydrogen_charges = [0.33818037500000003,
                    0.3322457708333334,0.34096977083333324,
                    0.32585266666666673,0.3194713958333334,0.3304236875,
                    0.31018468749999994,0.31812287499999997,0.31526810416666673,0.3455277708333333,
                    0.30178687499999995,0.31139518750000006,0.31443189583333336,0.3130669791666667,0.3360273541666668,
                    0.30178687499999995,0.31139518750000006,0.31443189583333336,0.3130669791666667,0.3360273541666668,0.34,
                    ]

hydrogen_charges = [-x + 1 for x in hydrogen_charges]

oxygen_charges = [7.326445125000002,
                  7.33217225,7.359637,
                  7.340767541666666,7.36089375,7.396371666666667,
                  7.369388583333333,7.3639823333333325,7.371821374999999,7.368475291666667,
                  7.3843305,7.377955208333333,7.371309041666668,7.3742492916666675,7.394824833333331,
                  7.3843305,7.377955208333333,7.371309041666668,7.3742492916666675,7.394824833333331,7.40,
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
