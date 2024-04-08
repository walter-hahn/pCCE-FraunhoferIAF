import numpy as np
import math
from scipy import sparse as sp

# Attempt to load custom settings from input.py if it exists
# Attempt to import max_size
try:
    from input import max_part_size as max_size
except ImportError:
    max_size = 4  # Default value

# Attempt to import min_size
try:
    from input import min_part_size as min_size
except ImportError:
    min_size = 4  # Default value

# Attempt to import no_systems
try:
    from input import number_systems as no_systems
except ImportError:
    no_systems = 1  # Default value

# Attempt to import number_interlaced
try:
    from input import number_internal_avg as number_interlaced
except ImportError:
    number_interlaced = 20  # Default value

# Attempt to import number_atoms
try:
    from input import number_spins as number_atoms
except ImportError:
    number_atoms = 144  # Default value

# Attempt to import base_concentrations
try:
    from input import base_concentrations, thickness, r_dipole_weight, t_max, usp_flag
except ImportError:
    base_concentrations = [1e-6]  # Default value
    thickness = 240  # Default value
    r_dipole_weight = 45  # Default value
    t_max = 500  # Default value
    usp_flag = 0

hbar = 6.626 * 10 ** -34 / (2 * math.pi)
mu0 = 1.257 * 10 ** -6
gamma_e = 28000
gamma_c = 10.705
Bz = 0.02
eZ = Bz * gamma_e
nZ = Bz * gamma_c
try:
    from input import gamma_b
except ImportError:
    gamma_b = gamma_e

Bx = 0.000001 * 50
Be = Bx * gamma_e
D = 2870
f = (D - eZ)
w = (D - eZ) * 2 * math.pi
hfs = 57
v0 = [0, 0, 1]
v1 = [1,1,1]/np.linalg.norm([1,1,1])
# spin matrices
sigma_z = sp.csr_matrix(np.array([[1, 0], [0, -1]]))
sigma_z_NV = np.array([[0, 0], [0, -2]])
sigma_x = sp.csr_matrix(np.array([[0, 1], [1, 0]]))
sigma_y = sp.csr_matrix(np.array([[0, -1j], [1j, 0]]))
sigma_zz = np.array([[1, 0], [0, 1]])

I_e = sp.csr_matrix(np.array([[1, 0], [0, 1]]))

# spin matrices spin 1
Sz = np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]])
Sx = 1 / np.sqrt(2) * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
Sy = 1j / np.sqrt(2) * np.array([[0, 1, 0], [-1, 0, 1], [0, -1, 0]])
Szz = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1]])

I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

S_plus = np.array([[0, 1], [0, 0]])
S_minus = np.array([[0, 0], [1, 0]])
S3 = []
S2 = []
S3.append(Sx)
S3.append(Sy)
S3.append(Sz)
S2.append(sigma_x)
S2.append(sigma_y)
S2.append(sigma_z)
# S2 = sp.csr_matrix(S2)

magx_base = S2[0] * 0.5
identity = sp.csr_matrix(np.identity(2))

x_gate = {}
gate = sigma_x
for i in range(12):
    gate = sp.kron(gate, I_e, format='csr')

    x_gate[i + 1] = sp.csr_matrix(gate)

exp_x = {}
exp_min_x = {}

for key, gate in x_gate.items():
    exp_x[key] = gate
    exp_min_x[key] = gate
