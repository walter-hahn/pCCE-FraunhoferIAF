import random
from constants import *
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse import linalg as spl


def get_hamiltonian(mf_positions, cluster_number, clusters, states, H_center, H_dict):
    """
    Generates the Hamiltonian for a given cluster within a spin system.

    This function constructs the Hamiltonian matrix for a specific cluster in a system of spins.
    The Hamiltonian includes interactions between spins in the cluster, interactions with a central
    spin, and bath-bath interactions among spins. The function accounts for dipole interactions and
    spin bath interactions using the provided matrices and positions.

    Parameters:
    mf_positions (list): List of magnetic field positions.
    cluster_number (int): The index of the current cluster in the `clusters` list.
    clusters (list): A list containing arrays of positions for each cluster.
    states (dict): A dictionary mapping positions to spin state matrices.
    H_center (ndarray): The Hamiltonian matrix of the central spin.
    H_dict (dict): A dictionary mapping positions to Hamiltonian matrices for defects.

    Returns:
    csr_matrix: The constructed sparse Hamiltonian matrix for the given cluster.

    Notes:
    The function uses Kronecker products to build up the Hamiltonian matrix, taking into account
    various interactions, including the central spin interactions and spin bath interactions.
    The Hamiltonian is computed as a sparse matrix for efficiency, especially useful for large systems.
    """
    positions = clusters[cluster_number]
    size = len(positions)
    p = [i for i in range(size)]
    # p = []
    # i = 0
    # for position in positions:
    #     p.append(i)
    #     i = i + 1
    element = tuple(p)
    size = len(element)

    # result = H_center
    result = csr_matrix(H_center)
    rel_positions = [positions[num] for num in element]

    # NV_Hamiltonian
    for s in range(size):
        result = sp.kron(result, I_e, format='csr')
    Cluster_H = csr_matrix((2 ** (size + 1), 2 ** (size + 1)))
    Cluster_H = Cluster_H + result

    # Logic for building the Hamiltonian is extracted from the original function
    # detect hamiltonians
    for d_number in element:
        result = I_e
        for number in element:
            if d_number == number:
                position = clusters[cluster_number][number]
                hamiltonian_defect = H_dict[str(position)]
                for m in range(len(mf_positions)):
                    if not np.any(np.all(mf_positions[m] == rel_positions, axis=1)):
                        dipole_interaction = compute_interaction_b(mf_positions[m] - position) #new compute interaction
                        matrix = states[str(mf_positions[m])]
                        if matrix[0, 0] == 1:
                            factor = 1
                        else:
                            factor = -1
                        hamiltonian_defect = hamiltonian_defect + sigma_z * dipole_interaction * factor / 4
                hamiltonian_defect = csr_matrix(hamiltonian_defect)
                result = sp.kron(result, hamiltonian_defect, format='csr')
            else:
                result = sp.kron(result, I_e, format='csr')
        Cluster_H = Cluster_H + result

    for d_number1 in element:
        for d_number2 in element:
            if d_number2 >= d_number1:
                if d_number1 == d_number2:
                    # interactions with center
                    dipole_interaction = compute_interaction_eb(positions[d_number2])
                    ham = []
                    # spin bath interaction
                    for i in range(1):

                        int_Hamiltonian = csr_matrix(np.array([1]))

                        for number in element:
                            if d_number2 == number:
                                int_Hamiltonian = sp.kron(int_Hamiltonian, S2[2], format='csr')

                            else:
                                int_Hamiltonian = sp.kron(int_Hamiltonian, I_e, format='csr')

                        int_Hamiltonian = sp.kron(sigma_z_NV, int_Hamiltonian, format='csr') * dipole_interaction / 4

                        ham.append(int_Hamiltonian)
                    result = np.zeros([ham[0].shape[1], ham[0].shape[1]])
                    for ele in ham:
                        result = result + ele
                    Cluster_H = Cluster_H + result

                else:
                    # bath-bath interactions
                    dipole_interaction = compute_interaction_b(positions[d_number2] - positions[d_number1])
                    ham = []

                    for i in range(3):
                        if i == 2:
                            factor = 4
                        else:
                            factor = -8

                        int_Hamiltonian = csr_matrix(np.array([1]))
                        for number in element:
                            if d_number1 == number:
                                int_Hamiltonian = sp.kron(int_Hamiltonian, S2[i], format='csr')

                            elif d_number2 == number:
                                int_Hamiltonian = sp.kron(int_Hamiltonian, S2[i], format='csr')

                            else:
                                int_Hamiltonian = sp.kron(int_Hamiltonian, I_e, format='csr')

                        int_Hamiltonian = sp.kron(I_e, int_Hamiltonian, format='csr') * dipole_interaction / factor

                        ham.append(int_Hamiltonian)

                    result = np.zeros([ham[0].shape[1], ham[0].shape[1]])
                    for ele in ham:
                        result = result + ele
                    Cluster_H = Cluster_H + result

    # Return the Hamiltonian
    Cluster_H = np.real(Cluster_H)
    Cluster_H = csr_matrix(Cluster_H)
    return Cluster_H


def compute_interaction_b(a):
    unit_vector = a / np.linalg.norm(a)
    dot_product = np.dot(v1, unit_vector)
    dipole_interaction = mu0 * (gamma_b * 10 ** 6) ** 2 * math.pi / (np.linalg.norm(a) * 10 ** -9) ** 3 * \
                         (1 - 3 * dot_product ** 2) * 10 ** -6 * hbar
    return dipole_interaction


def compute_interaction_eb(a):
    unit_vector = a / np.linalg.norm(a)
    dot_product = np.dot(v1, unit_vector)
    dipole_interaction = mu0 * (gamma_e * 10 ** 6) * (gamma_b * 10 ** 6) * math.pi / (np.linalg.norm(a) * 10 ** -9) ** 3 * \
                        (1 - 3 * dot_product ** 2) * 10 ** -6 * hbar
    return dipole_interaction


def get_time_propagator(Cluster_H, tau):
    # Convert Cluster_H to CSC format for efficiency
    Cluster_H_csc = csc_matrix(Cluster_H)

    # Compute the matrix exponential
    matrix_exp = spl.expm(-1j * Cluster_H_csc * tau)

    # Convert the result back to CSR format if needed
    matrix_exp_csr = csr_matrix(matrix_exp)

    return matrix_exp_csr


def time_propagation(matrix_exp, psi0):
    """
        Propagates a quantum state in time using a given matrix exponential.

        Parameters:
        matrix_exp (ndarray): The matrix exponential of the Hamiltonian representing time evolution.
        psi0 (ndarray): The initial quantum state vector.

        Returns:
        ndarray: The quantum state after propagation (psi_tau).
    """
    psi_tau = matrix_exp.dot(np.transpose(psi0))
    return psi_tau


def get_probabilities(matrix_exp_step, positions, magx, num_states, testing=0):
    """
    Calculates the average probability of quantum state evolution under a given matrix.

    This function computes the average probability of a quantum state after undergoing
    evolution as dictated by the provided matrix exponential step. It generates random
    initial states and calculates the probability of these states evolving to a specific
    state under the Hamiltonian dynamics.

    Parameters:
    matrix_exp_step (ndarray): The matrix exponential step for state evolution.
    positions (list): The positions of spins or particles in the system.
    magx (ndarray): The magnetic field matrix.
    num_states (int): Number of random states to generate for averaging.
    testing (int, optional): A flag for testing purposes (default is 0).
                             If set, uses a fixed seed for random number generation.

    Returns:
    float: The average probability computed over `num_states` random initial states.
    """
    len_positions = len(positions)
    state_size = 2 ** (len_positions + 1)

    # Precompute outside the loop
    exp_x_positions = exp_x[len_positions]
    if testing:
        np.random.seed(42)

    ran1 = np.random.uniform(0, 1, size=(num_states, 2 ** len_positions))
    ran2 = np.random.uniform(0, 1, size=(num_states, 2 ** len_positions))

    w1 = np.sqrt(-2 * np.log(ran1)) * np.cos(2 * ran2 * np.pi)
    w2 = np.sqrt(-2 * np.log(ran1)) * np.sin(2 * ran2 * np.pi)

    factors = w1 + 1j * w2

    norm = 2 * np.sum(w1 ** 2 + w2 ** 2, axis=1, keepdims=True)

    state1 = np.zeros((num_states, state_size), dtype=complex)
    state1[:, :2 ** len_positions] = factors
    state1[:, 2 ** len_positions:] = factors

    psi0 = state1 / np.sqrt(norm)

    psi_tau = matrix_exp_step @ psi0.T
    psi_x = exp_x_positions @ psi_tau
    psi_two_tau = matrix_exp_step @ psi_x

    prob1 = np.real(np.einsum('ij,ji->i', psi_two_tau.T.conj(), magx @ psi_two_tau)) #check values on this line
    prob2 = 0.5

    prob = prob1 / prob2

    return np.mean(prob)


def solve_hahn(mf_positions, t_max, n_steps, clusters, states, H_center, H_dict, test=0):
    """
        Solves the Hahn echo problem for a system of spins.

        This function computes the probabilities of quantum states over time for a spin system,
        particularly for the Hahn echo experiment. It involves generating the Hamiltonian for
        each cluster, propagating the states in time, and computing probabilities.

        Parameters:
        mf_positions (list): List of magnetic field positions.
        t_max (float): The maximum time for the Hahn echo experiment.
        n_steps (int): The number of time steps for the experiment.
        clusters (list): A list of clusters, each containing positions of spins.
        states (dict): A dictionary mapping positions to spin state matrices.
        H_center (ndarray): The Hamiltonian of the central spin.
        H_dict (dict): A dictionary of Hamiltonians for different positions.
        test (int, optional): A flag for testing purposes (default is 0).

        Returns:
        tuple: A tuple containing two elements: `time_prob` and `all_probs`.
               `time_prob` is an array of time points, and `all_probs` is an array of corresponding probabilities.
    """
    if test == 1:
        random.seed(48)
    all_probs = np.zeros([n_steps - 1, len(clusters)])
    t_step = t_max / n_steps
    tau = t_step / 2

    mag_x = magx_base.copy()
    mag_2x = magx_base.copy()

    for _ in range(max_size):
        mag_x = sp.kron(mag_x, identity, format='csr')
    for _ in range(max_size * 2):
        mag_2x = sp.kron(mag_2x, identity, format='csr')

    for cluster_number in range(len(clusters)):
        Cluster_H = get_hamiltonian(mf_positions, cluster_number, clusters, states, H_center, H_dict)
        # start_time_tp = time.time_ns()
        matrix_exp = get_time_propagator(Cluster_H, tau)
        # finish_time_tp = time.time_ns()
        # print(finish_time_tp-start_time_tp)
        matrix_exp_step = matrix_exp.copy()

        positions = clusters[cluster_number]

        if len(positions) == max_size:
            magx = mag_x
        elif len(positions) == (max_size*2):
            magx = mag_2x
        else:
            # print("generating new magx")
            magx = magx_base
            for _ in range(len(positions)):
                magx = sp.kron(magx, identity, format='csr')

        current_step = -1
        num_states = 10 if Cluster_H.shape[1] > (2 ** (1 + max_size) + 2) else 32
        for t in np.linspace(t_step * 2, t_max, n_steps - 1):
            current_step += 1
            matrix_exp_step = matrix_exp_step.dot(matrix_exp)
            prob = get_probabilities(matrix_exp_step, positions, magx, num_states)
            all_probs[current_step, cluster_number] = prob
    t = t_step
    time_prob = [[t + i * t_step, np.prod(element)] for i, element in enumerate(all_probs)]
    time_prob = np.stack(time_prob, axis=0)
    return time_prob, all_probs

