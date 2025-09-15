import matplotlib.pyplot as plt
import matplotlib_parameters

# Data
k_grid = ['1x1x1', '2x2x1', '3x3x1', '4x4x1', '5x5x1']
wf_values = [4.8853, 5.0609, 5.1418, 5.1418, 5.1531]
experimental_wf = 5.33

# Plotting
fig, ax = plt.subplots()
bars = ax.bar(k_grid, wf_values, color='skyblue', edgecolor='black')

# Horizontal experimental line
ax.axhline(experimental_wf, color='darkgreen', linestyle='--', linewidth=2)
ax.text(len(k_grid) - 1, experimental_wf + 0.02, 'Exp.', color='darkgreen', fontsize=12)

# Labels and formatting
ax.set_ylabel('WF (eV)')
ax.set_xlabel('k-point grid')
ax.set_ylim(4.8, 5.4)

plt.show()

