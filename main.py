import time
import numpy as np
import argparse
from mpi4py import MPI
from utils import build_lattice, constrained_clustering, mf_bath, get_time_prob
from solve_hahn import solve_hahn
from data_manager import save_result_data, save_avg_result_data, load_data
from constants import *
import sys

# Attempt to override constants with values from input.py, if exists
try:
    from input import *
except ImportError:
    pass  # If input.py doesn't exist, just use constants from constants.py


# Define command-line arguments
def get_command_line_args():
    parser = argparse.ArgumentParser(description="Simulation Parameters")
    # Define arguments for each configurable parameter, matching those in input.py and constants.py
    parser.add_argument("--info", action="store_true", help="Display information about input methods")
    parser.add_argument("--max_size", type=int, help="Maximum cluster size", default=max_size)
    parser.add_argument("--min_size", type=int, help="Minimum cluster size", default=min_size)
    parser.add_argument("--no_systems", type=int, help="Number of systems", default=no_systems)
    parser.add_argument("--t_max", type=int, help="Maximum time for simulation", default=t_max)
    parser.add_argument("--number_atoms", type=int, help="Number of atoms in a system", default=number_atoms)
    parser.add_argument("--number_interlaced", type=int, help="Number of interlaced simulations",
                        default=number_interlaced)
    parser.add_argument("--base_concentrations", nargs='+', type=float, help="Base concentrations", default=[1e-6])
    parser.add_argument("--thickness", type=int, help="Layer thickness of the spin distribution", default=240)
    parser.add_argument("--r_dipole_weight", type=int, help="Dipole weight", default=45)
    parser.add_argument("--gamma_b", type=int, help="Gamma B", default=28000)
    parser.add_argument("--usp_flag", type=int, help="USP Flag", default=0)
    return parser.parse_args()


args = get_command_line_args()

# Override values with command-line arguments This is redundant if defaults are set in the parser, as argparse will
# use the default value if not specified by the user
if args.info:
    print("There are two ways to provide inputs for the simulation:\n"
          "1. Command-Line Arguments: Directly pass parameters when running the script.\n"
          "   Example: python main.py --no_systems 10 --number_atoms 100\n"
          "2. Modifying input.py: Create or modify input.py with your desired parameters.\n"
          "   This file should define variables like no_systems, number_atoms, etc.")
    sys.exit(0)
no_systems = args.no_systems
number_atoms = args.number_atoms
t_max = args.t_max
number_interlaced = args.number_interlaced
max_size = args.max_size
min_size = args.min_size
base_concentrations = args.base_concentrations
thickness = args.thickness
r_dipole_weight = args.r_dipole_weight
gamma_b = args.gamma_b
usp_flag = args.usp_flag

# Ensure this script runs only as the main program
if __name__ == '__main__':
    # Initialize MPI communication
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    # Start timing the execution
    start_main = time.time()
    print("Starting simulation....")

    # Prepare parameter list for simulations
    parameter_list = [(conc, thickness) for conc in base_concentrations]

    # Iterate over each set of parameters
    for concentration, thickness in parameter_list:
        # Calculate dipole radius based on concentration and user defined weights
        r_dipole = r_dipole_weight * ((1 * 10 ** -6) / concentration) ** (1 / 3) * 3 ** (1 / 3)
        # Adjust r_dipole based on thickness
        if thickness < (2 * r_dipole):
            r_dipole *= (2 * r_dipole / thickness) ** 0.5
        results = []

        # Adjust time factor based on thickness
        t_factor = 150 / thickness if thickness < 150 else 1

        # Debugging print statement
        print(size, rank, concentration, thickness)

        # Simulate for the number of systems
        for i in range(no_systems):
            # systems = [load_data("pickleFiles/positions_20_320230120.pkl")]
            # mf_positions = load_data("pickleFiles/positions_20_320230120.pkl")
            systems, mf_positions = build_lattice(concentration, number_atoms, thickness)
            final_time_prob = []

            for all_positions in systems:
                # Determine number of clusters based on positions and max cluster size
                n_clusters = len(all_positions) // max_size
                if len(all_positions) < 30:
                    r_dipole = 10000
                clusters, num_large_clusters = constrained_clustering(all_positions, n_clusters, r_dipole)
                # clusters, num_large_clusters = constrained_clustering_agg(all_positions, n_clusters, r_dipole)
                all_all_probs = np.zeros_like(all_positions)  # Initialize probability array

                for j in range(number_interlaced):

                    # Generate magnetic field and Hamiltonian for the system
                    states, H_center, H_dict = mf_bath(all_positions, mf_positions)  # system is always 0
                    # Solve the Hahn echo problem
                    time_prob, all_probs = solve_hahn(mf_positions, t_max, 40, clusters, states, H_center, H_dict)

                    # Accumulate probabilities
                    if j == 0:
                        all_all_probs = all_probs
                        probs_large = all_all_probs
                    else:
                        all_all_probs += all_probs

                final_all_probs = all_all_probs / number_interlaced  # Average probabilities
                single_prob = get_time_prob(final_all_probs, num_large_clusters, t_max)  # Extract single probability
                final_time_prob.append(single_prob)  # Collect probabilities for all systems
            results.append(final_time_prob.copy())  # Collect results for all concentrations

        rslts_list = comm.gather(results, root=0)  # Gather results across MPI nodes
        mean_result = np.mean(results, axis=0)  # Compute mean result for the current set of parameters

        if rank == 0:
            # Save results and mean results to file
            save_result_data(rslts_list, concentration, thickness)
            save_avg_result_data(mean_result)

    # End timing and print elapsed time
    end_main = time.time()
    print("Time elapsed:", end_main - start_main)
