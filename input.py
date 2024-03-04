# This file allows for custom configuration of the simulation parameters.
# Modify these values as needed for your simulation.

max_part_size = 4    # max_part_size
min_part_size = 4    # min_part_size

number_systems = 1  # number_systems
t_max = 500
number_internal_avg = 20  # number_internal_avg
number_spins = 144      # number_spins
base_concentrations = [1e-6]
thickness = 240

r_dipole_weight = 45
# Used to modify r_dipole
# r_dipole = r_dipole_weight * ((1 * 10 ** -6) / concentration) ** (1 / 3) * 3 ** (1 / 3)

