import time
import argparse
from mpi4py import MPI
from utils import build_lattice, constrained_clustering, mf_bath, get_time_prob
from solve_hahn import solve_hahn
from data_manager import save_result_data, save_avg_result_data, load_data
from constants import *
import matplotlib.pyplot as plt
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
    parser.add_argument("--base_concentrations", nargs='+', type=float, help="Base concentrations", default=base_concentrations)
    parser.add_argument("--thickness", type=int, help="Layer thickness of the spin distribution", default=thickness)
    parser.add_argument("--r_dipole", type=int, help="Dipole weight", default=r_dipole)
    parser.add_argument("--time_step", type=int, help="Time step(ms)", default=time_step)
    parser.add_argument("--gamma_b", type=int, help="Gamma B", default=gamma_b)
    parser.add_argument("--hf_for_P1", type=int, help="Flag to change spin bath from Electron spins to P1 centres", default=hf_for_P1)
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
r_dipole = args.r_dipole
gamma_b = args.gamma_b
hf_for_P1 = args.hf_for_P1

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
        results = []
        # Adjust time factor based on thickness
        t_factor = 150 / thickness if thickness < 150 else 1

        # Debugging print statement
        print(size, rank, concentration, thickness)

        # Simulate for the number of systems
        for i in range(no_systems):
            systems, mf_positions = build_lattice(concentration, number_atoms, thickness, r_dipole)
            final_time_prob = []

            for all_positions in systems:
                # Determine number of clusters based on positions and max cluster size
                n_partitions = int(len(all_positions) // ((max_size+min_size)/2))
                if hf_for_P1:
                    if len(all_positions) < 30:
                        r_dipole = 10000
                clusters, num_large_clusters = constrained_clustering(all_positions, n_partitions, r_dipole)
                all_all_probs = np.zeros_like(all_positions)  # Initialize probability array

                for j in range(number_interlaced):

                    # Generate magnetic field and Hamiltonian for the system
                    states, H_center, H_dict = mf_bath(all_positions, mf_positions)  # system is always 0
                    # Solve the Hahn echo problem
                    time_prob, all_probs = solve_hahn(mf_positions, t_max, time_step, clusters, states, H_center, H_dict)

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

    # End timing and print elapsed time
    end_main = time.time()
    print("Time elapsed:", end_main - start_main)
    if rank == 0:
        P1_list_1 = []
        for element in rslts_list:
            for elem in element:
                combined_result = []
                if hf_for_P1:
                    # Compute product of results from multiple systems
                    sys_list = [elem[i] for i in range(5)]  # Get systems from 0 to 4
                    for i in range(len(sys_list[0])):
                        result = np.prod([sys[i, 1] for sys in sys_list])  # Using numpy product for elegance and speed
                        combined_result.append([sys_list[0][i, 0], result])
                else:
                    # Use only the first system's results
                    sys0 = elem[0]
                    for i in range(len(sys0)):
                        result = sys0[i, 1]
                        combined_result.append([sys0[i, 0], result])
                P1_list_1.append(combined_result)
        # Compute the average of the results across all entries in P1_list_1
        mean_result = np.average(P1_list_1, axis=0)
    save_avg_result_data(mean_result)
    save_result_data(P1_list_1, concentration, thickness)
    plt.figure(figsize=(10, 5))
    plt.plot(mean_result[:, 0], mean_result[:, 1], marker='o', linestyle='-')
    plt.title('Plot of Mean Result')
    plt.savefig('mean_result_plot.png')
    plt.xlabel('time $2\\tau$ (μs)')  # using LaTeX for subscript and special characters
    plt.ylabel('$M_x$')  # using LaTeX for subscript
    plt.show()
