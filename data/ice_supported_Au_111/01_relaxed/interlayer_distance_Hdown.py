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

int_dist = [0.07699676287325538,
            0.07666166886860934,0.07654223002018967,
            0.07585741202236826,0.07617577733187392,0.07577744627055744,
            0.07569779979964236,0.07598399276901163,0.07592492717611182,0.07540288473187934,
            0.07575897550474237,0.07605681236174339,0.07601045046390392,0.0760245660463787,0.07545964761111257,
            ]


int_dist = [x * z_cell for x in int_dist]

fig, ax = plt.subplots()
ax.scatter(number_bilayer, int_dist, color='blue', s=100)
ax.set_xlabel("Number of bilayers")
ax.set_ylabel("Interlayer distance ($\mathrm{\AA}$)")
#ax.set_xlim([-15,5])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(1))
plt.show()

