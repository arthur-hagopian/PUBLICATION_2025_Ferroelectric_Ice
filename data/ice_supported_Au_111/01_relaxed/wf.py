#!/usr/bin/python3.10
# Arthur Hagopian <arthur.hagopian@umontpellier.fr>, version 01/10/2024

import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

number_bilayers = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]

wf = [1.79,1.80,1.80,2.29,2.17,1.39,4.88,6.18,8.34,8.36,8.37,8.38,8.38]

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

