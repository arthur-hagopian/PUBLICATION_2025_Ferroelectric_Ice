#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

z_cell = 49.2004013061999999

number_bilayer = [2,
                  3,3,
                  4,4,4,
                  5,5,5,5,
                  6,6,6,6,6
                  ]

int_dist = [0.08157601354581973,
            0.07995860472403873,0.07871310055166647,
            0.07829666817053815,0.07849391916851167,0.07767545994445313,
            0.07692249272760301,0.07622574404753374,0.07623295039464084,0.07640642837222156,
            0.07660066505847365,0.07606955339481158,0.07611261383471696,0.07625915638816055,0.07636296909051832
            ]

colors = ['blue',
          'blue','red',
          'blue','red','green',
          'blue','red','green','orange',
          'blue','red','green','orange','purple',
          'blue','red','green','orange','purple','pink',
          ]

int_dist = [x * z_cell for x in int_dist]

fig, ax = plt.subplots()
i = 0
while i < len(number_bilayer):
    ax.scatter(number_bilayer[i], int_dist[i], color=colors[i], s=100)
    i += 1
ax.set_xlabel("Number of bilayers")
ax.set_ylabel("Interbilayer distance ($\mathrm{\AA}$)")
#ax.set_xlim([-15,5])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(1))
plt.show()

