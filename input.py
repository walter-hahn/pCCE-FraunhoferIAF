# This file allows for custom configuration of the simulation parameters.
# Modify these values as needed for your simulation.

max_part_size = 4                     # Maximum partition size
min_part_size = 4                     # Minimum partition size

number_systems = 1                    # Number of systems considered. Each system corresponds to a random spatial distribution of spins. MPI parallelization distributes theses systems among available cores.
t_max = 500                           # Maximum time in µs
number_internal_avg = 1               # Number of averages performed internally (Monte Carlo bath state sampling)
number_spins = 144                     # Total number of bath-spins included in the calculation
concentration = 1e-6                  # Concentration (Fraction of atoms replaced by defects)
thickness = 240                       # Diamond-layer thickness along the z-direction in nm. The central spin is positioned at the center of the layer.
hf_for_P1 = 0                         # Set usp_flag=1 to simulate P1-centers, usp_flag=0 to simulate electron spins in the bath.
gamma_b = 28000                       # Gyromagnetic ratio of the bath spins in MHz/T. Default is for electron spins.
r_dipole = 45                         # Dipole-radius of the bathspins defining the distance in which dipolar interactions between bath spins are considered.
time_step = 10                        # Number of time steps to t_max