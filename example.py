# Example
#
# Consider a 24-ply quasi-isotropic laminate. The ply properties and ply
# thickness are listed in this file. This example file will illustrate how
# to use the functions in this repo to determine:
# - Engineering constants of the laminate
# - Deformation and stress distribution due to an applied load
# - Deformation and stress distribution because of a temperature change

import clt
import numpy as np

# Layup
t = 0.15E-3
n = 24
layup = clt.QI_layup(n)
z = clt.ply_edges(t, n)

# material properteis
E1, E2, v12, G12 = 120E9, 10E9, 0.28, 5E9
C = clt.stiffness_matrix(E1, E2, v12, G12)
alpha = np.array([0.2E-6, 30E-6, 0])  # CTE in material CS

# ABD and abd matrices
C_r = clt.rotate_C(C, layup)  # List with ply stiffness matrices
ABD = clt.ABD_matrix(C_r, z)
abd = np.linalg.inv(ABD)

# Laminate engineering constants
E_x = 1/(abd[0, 0]*n*t)
E_y = 1/(abd[1, 1]*n*t)
G_xy = 1/(abd[2, 2]*n*t)

# Deformation and stresses due to a load applied in x direction
F = np.array([1E5, 0, 0, 0, 0, 0])
d = abd@F
stress_r, z_int = clt.ply_stress(d, C_r, z)
stress = clt.rotate_stress_to_matCS(stress_r, layup)

# Plot results
clt.plot_stress(stress_r, z_int, 0)  # plot stress in the 1* direction

# Deformation and stresses due to cooling from 220 to room temperature
deltaT = -200
alpha_r = clt.rotate_alpha(alpha, layup)
Fth = clt.thermal_force(C_r, alpha_r, z, deltaT)
d = abd@Fth
stress_r, z_int = clt.ply_stress(d, C_r, z, alpha_r, deltaT)
stress = clt.rotate_stress_to_matCS(stress_r, layup)

# Plot results
clt.plot_stress(stress, z_int, 0)  # plot stress in the 1 direction
