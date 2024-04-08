# This file allows for custom configuration of the simulation parameters.
# Modify these values as needed for your simulation.

max_part_size = 1
min_part_size = 1
number_systems = 2
t_max = 500
number_internal_avg = 1
number_spins = 20
base_concentrations = [1e-6]
thickness = 240         # Layer thickness of the spin distribution


usp_flag = 1            # Change spin bath from Electron spins to P1 centres
gamma_b = 28000         # used only when usp_flag is set to 1
r_dipole_weight = 45
# Used to modify r_dipole
# r_dipole = r_dipole_weight * ((1 * 10 ** -6) / concentration) ** (1 / 3) * 3 ** (1 / 3)

