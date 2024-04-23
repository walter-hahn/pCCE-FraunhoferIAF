import random
from k_means_constrained import KMeansConstrained
from constants import *


def generate_position(iterations_x, iterations_y, iterations_z, var1, var2, var3, choice_list):
    """
        Generates a random position vector within a lattice structure.

        This function creates a random position within a defined lattice space.
        It combines base lattice vectors with random coefficients to generate the position.

        Parameters:
        iterations_x, iterations_y, iterations_z (int): The range of iterations in each dimension.
        var1, var2, var3 (array_like): Base vectors for the lattice in x, y, and z dimensions.
        choice_list (list): A list of initial base positions to choose from.

        Returns:
        array_like: A randomly generated position vector within the lattice.
    """
    rand1 = random.randint(-iterations_x, iterations_x)
    rand2 = random.randint(-iterations_y, iterations_y)
    rand3 = random.randint(-iterations_z, iterations_z)
    position_vector = random.choice(choice_list) + rand1 * var1 + rand2 * var2 + rand3 * var3
    return position_vector / 40 * 3.567


def build_lattice(concentration, number_spins, thickness, r_dipole):
    """
        Constructs a lattice of spins based on specified parameters.

        This function builds a lattice structure with a given concentration of spins,
        considering a specific number of spins and the thickness of the lattice layer.
        It also considers the radius of interaction for bath spins.

        Parameters:
        concentration (float): The concentration of spins in the lattice.
        number_spins (int): The total number of spins in the lattice.
        thickness (float): The thickness of the lattice layer.

        Returns:
        tuple: A tuple containing the arrangement of spins in the lattice and the positions
               of mean field (mf) spins.
    """
    b1 = np.array([0, 0, 0])
    b2 = np.array([0, 2, 2])
    b3 = np.array([2, 0, 2])
    b4 = np.array([2, 2, 0])
    b5 = np.array([3, 3, 3])
    b6 = np.array([3, 1, 1])
    b7 = np.array([1, 3, 1])
    b8 = np.array([1, 1, 3])

    var1 = np.array([4, 0, 0])
    var2 = np.array([0, 4, 0])
    var3 = np.array([0, 0, 4])
    base_list = (b1, b2, b3, b4, b5, b6, b7, b8)

    if thickness < (2 * r_dipole):
        r_dipole = r_dipole * (2 * r_dipole / thickness) ** (1 / 2)
        print("Warning!! r_dipole as been modified w.r.t the thickness")
        print("r_dipole = " + str(r_dipole))

    # create random lattice sites

    size_factor = (10 ** -6 / concentration) ** (1 / 3) 
    if thickness < 120*(10 ** -6 / concentration) ** (1 / 3):
        size_factor = size_factor * (120*(10 ** -6 / concentration) ** (1 / 3) / thickness) ** (1 / 2)
    iterations_x = int(500 * size_factor)
    iterations_y = int(500 * size_factor)
    iterations_z = int(500 * size_factor)

    number_atoms = 8 * 8 * iterations_x * iterations_y * iterations_z
    number_N = int(concentration * number_atoms)

    all_positions = []
    for _ in range(number_N):
        real_position = generate_position(iterations_x, iterations_y, iterations_z, var1, var2, var3, base_list)
        while any(np.array_equal(element[:3], real_position[:3]) for element in all_positions) or np.array_equal(
                real_position[:3], [0, 0, 0]):
            real_position = generate_position(iterations_x, iterations_y, iterations_z, var1, var2, var3, base_list)
        if abs(np.dot(real_position, v1)) < thickness / 2:
            all_positions.append(real_position)

    all_posit = sorted(all_positions, key=np.linalg.norm)
    all_pos = all_posit[:number_spins]
    last_atom = len(all_pos) + 1
    if hf_for_P1:
        choices = [0, 1, 2, 3, 4]
        decisions = random.choices(choices, weights=(1 / 3, 1 / 4, 1 / 4, 1 / 12, 1 / 12), k=len(all_pos))
        assignment = [[] for _ in range(5)]
    else:
        choices = [0]
        decisions = random.choices(choices, k=len(all_pos))
        assignment = [[] for _ in range(1)]

    for i, decision in enumerate(decisions):
        assignment[decision].append(all_pos[i])

    for assign_list in assignment:
        missing_atoms = max_part_size - len(assign_list) % max_part_size
        for _ in range(missing_atoms):
            assign_list.append(all_posit[last_atom])
            last_atom += 1

    assignment = [np.stack(assign_list) for assign_list in assignment]

    max_dist = np.linalg.norm(all_pos[-1])
    mf_positions = [pos for pos in all_posit if np.linalg.norm(pos) < max_dist + r_dipole * (2 / 3)]
    mf_positions = np.stack(mf_positions, axis=0)

    return assignment, mf_positions


def get_constrained_clusters(positions, n_partitions):
    """
        Groups positions into a specified number of clusters with constraints.

        This function applies KMeans clustering with constraints on the size of each cluster.
        It groups the given positions into a specified number of clusters.

        Parameters:
        positions (array_like): An array of positions to be clustered.
        n_partitions (int): The number of clusters to form.

        Returns:
        tuple: A tuple containing the labels for each position and the cluster centers.
        """
    clf = KMeansConstrained(
        n_clusters=n_partitions,
        size_min=min_part_size,
        size_max=max_part_size,
        random_state=0
    )
    clf.fit_predict(positions)

    return clf.labels_, clf.cluster_centers_


def constrained_clustering(all_positions, n_partitions, r_dipole):
    a, centers = get_constrained_clusters(all_positions, n_partitions)
    print(a)
    # sort position vectors
    temp = []
    exact_clusters = {}

    for i in range(len(centers)):
        for j in range(len(all_positions)):

            if a[j] == i:
                temp.append(all_positions[j])
        temp = np.stack(temp, axis=0)
        exact_clusters[i] = temp
        temp = []

    joined_clusters = []
    half_clusters = []

    num_multiply_clusters = 0

    subclusters = list(exact_clusters.values())
    for i in range(len(subclusters)):
        for j in range(len(subclusters)):
            if j > i:
                if i != j:
                    if np.linalg.norm(centers[i] - centers[j]) < r_dipole:
                        # Merge clusters if they are closer than the dipole radius
                        cluster_1 = subclusters[i]
                        cluster_2 = subclusters[j]
                        con = np.concatenate((cluster_1, cluster_2))
                        joined_clusters.append(con)
                        half_clusters.append(cluster_1)
                        half_clusters.append(cluster_2)
    clusters = {}
    i = 0

    for joined_cluster in joined_clusters:
        clusters[i] = joined_cluster
        i = i + 1
        num_multiply_clusters = num_multiply_clusters + 1

    for subcluster in subclusters:
        clusters[i] = subcluster
        i = i + 1
        num_multiply_clusters = num_multiply_clusters + 1

    for half_cluster in half_clusters:
        clusters[i] = half_cluster
        i = i + 1

    return clusters, num_multiply_clusters


def mf_bath(all_positions, mf_positions):
    """
        Generates a magnetic field bath based on given positions.

        This function creates a bath of magnetic fields for a set of positions,
        defining the Hamiltonians for the center and defects in the bath.

        Parameters:
        all_positions (array_like): An array of all positions in the system.
        mf_positions (array_like): Positions where magnetic fields are applied.

        Returns:
        tuple: A tuple containing the states, central Hamiltonian, and Hamiltonian dictionary.
        """

    # Hamiltonian
    H_center = np.zeros((3, 3))  # eZ*Sz+D*Szz#+Be*Sx/(2*np.sqrt(2))
    H_center = H_center[1:3, 1:3]

    H_defect = np.zeros((2, 2))  # (eZ)*sigma_z*1/2

    # Create states dictionary
    num_positions = len(mf_positions)
    random_indices = np.random.randint(2, size=num_positions)
    state_list = [np.array([[1, 0], [0, 0]]), np.array([[0, 0], [0, 1]])]
    states = {str(po): state_list[idx] for po, idx in zip(mf_positions, random_indices)}

    # Create H_dict
    H_dict = {str(pos): H_defect for pos in all_positions}  # +random.choice(level_shift)*sigma_z*1/2

    return states, H_center, H_dict


def get_time_prob(all_probs, num_large_clusters, t_max):
    """
        Calculates time-dependent probabilities for a system with multiple clusters.

        This function computes probabilities over time for a system composed of
        various clusters, plotting the results as a function of time.

        Parameters:
        all_probs (array_like): An array of probabilities for different clusters.
        num_large_clusters (int): The number of large clusters in the system.
        t_max (float): The maximum time for which probabilities are calculated.

        Returns:
        array_like: A 2D array with time and corresponding probabilities.
        """
    t_range = np.linspace((2 * t_max / time_step), t_max, time_step - 1)
    probabilities = []

    # Vectorized approach to calculate probabilities
    for element in all_probs:
        prod_first_elements = np.prod(element[:num_large_clusters])
        if len(element) > num_large_clusters:
            prod_other_elements = np.prod(1.0 / np.array(element[num_large_clusters:]))
        else:
            prod_other_elements = 1
        probabilities.append(prod_first_elements * prod_other_elements)

    time_prob = np.column_stack((t_range, probabilities))

    return time_prob
