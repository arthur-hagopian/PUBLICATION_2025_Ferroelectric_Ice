# Re-import necessary libraries since execution state was reset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib_parameters

# Define the new range of values for A and p
A_values = np.linspace(9, 20, 100)  # 100 points between 9 and 150
p_values = np.linspace(1, 2, 100)  # 100 points between 1 and 2

# Create a new meshgrid with the updated A and p ranges
A, p = np.meshgrid(A_values, p_values)

# Compute Delta Phi using the given relation
Delta_Phi = (p * 3.33564e-30) / ((A * 1e-20) * 8.854e-12)

# Plot the 2D colored map with the last used colormap (coolwarm)
plt.figure(figsize=(8, 6))
contour = plt.contourf(A, p, Delta_Phi, levels=100, cmap='coolwarm')
cbar = plt.colorbar(contour)
cbar.set_label(r'$\Delta \Phi_z$ (V)')

# Labels and title
plt.xlabel(r'Molecular occupancy area ($\AA$)')
plt.ylabel(r'Vertical dipole p$_z$ (D)')

# Show the plot
plt.show()

