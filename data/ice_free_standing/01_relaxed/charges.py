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

hydrogen_charges = [0.3495789166666667,
                    0.2997377916666666, 0.3191810833333333,
                    0.32414189583333336, 0.3054688125000001, 0.32558606250000005,
                    0.3369928333333334, 0.3202915624999999, 0.6036403125, 0.34092410416666663,
                    0.48977658333333346, 0.7019136041666667, 0.34520027083333327, 0.38494681249999996, 0.40717129166666677,
                    0.31259725000000016, 0.3000140833333333, 0.28297947916666677, 0.3056162083333333, 0.3978047291666667, 0.5480007708333335,
                    ]

hydrogen_charges = [-x + 1 for x in hydrogen_charges]

oxygen_charges = [7.300844958333332,
                  7.394143791666667, 7.3680205416666675,
                  7.378535916666667, 7.386291458333333, 7.3247799166666665,
                  7.392441375, 7.389998166666666, 6.661443666666667, 7.352413166666667,
                  6.984060333333333, 6.681616208333332, 7.264963333333333, 7.2788208333333335, 7.132521999999999,
                  7.502935083333334, 7.397097916666669, 7.369020625, 7.367045083333334, 7.102133833333334, 6.967748,
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
axs[2].set_xlabel('Number of Bilayers')
axs[2].set_ylabel('Charge (e)')
axs[2].xaxis.set_minor_locator(AutoMinorLocator())

for i in range(len(number_bilayer)):
    marker_style = 'D' if i % 2 == 0 else '*'  # Alternate between diamond and star markers
    axs[2].scatter(number_bilayer[i], hyd_plus_ox[i], color=colors[i], s=100)

# Show plot
plt.show()
