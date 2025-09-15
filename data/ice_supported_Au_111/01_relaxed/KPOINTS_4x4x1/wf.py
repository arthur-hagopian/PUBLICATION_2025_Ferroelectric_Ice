#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

number_bilayers = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]

wf = [2.22,2.2733,2.20,2.1894,1.9729,1.5993,5.1418,6.4087,8.3499,8.386800000000001,8.3753,8.3741,8.3797]

# Plot
fig, ax =plt.subplots()
i = 0
x = number_bilayers
y = wf
ax.scatter(x, y, s=200)

ax.set_xlabel('Number of bilayers (H-up / H-down)')
ax.set_ylabel('WF (eV)')

ax.set_xlim([-7,7])
#ax.set_ylim(bottom=0)
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(2))

ax.legend()
plt.show()

