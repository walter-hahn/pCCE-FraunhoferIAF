
# partition-based Cluster-Correlation Expansion (pCCE) method

This is an implementation of the pCCE method proposed in the manuscript arXiv:2401.16169 (https://arxiv.org/abs/2401.16169). The main focus is on the calculation of the spin-coherence decay of the NV-center spin in a bath of P1-defect spins. The code can also be used in a broader context.


## Features

- Creation of Hamiltonian matrices for spin clusters
- Time propagation of quantum states
- Calculation of probabilities for Hahn echo experiments
- Customizable simulation parameters via `constants.py`
- Parallelization using MPI

## Installation

To install and run the project:

1. **Clone the Repository:**
git clone https://github.com/walter-hahn/pCCE-FraunhoferIAF cd pCCE

2. **Install Dependencies:**
pip install -r requirements.txt

## Usage

Execute the main script to run the simulation:

```bash
python main.py
```
### Configuration
There are two ways to configure simulation parameters:
1. **Modifying input.py**: For persistent changes, you can modify the input.py file. This file overrides the default settings specified in constants.py. If input.py does not exist, the simulation will fall back to using the default values.

Sample input.py configuration:
```bash
max_part_size = 4                     # Maximum partition size    
min_part_size = 4                     # Minimum partition size
number_systems = 1                    # Number of systems considered. Each system corresponds to a random spatial distribution of spins. MPI parallelization distributes theses systems among available cores.
t_max = 500                           # Maximum time in µs
number_internal_avg = 20              # Number of averages performed internally (Monte Carlo bath state sampling)
number_spins = 144                    # Total number of bath-spins included in the calculation       
concentration = 1e-6                  # Concentration (Fraction of atoms replaced by defects)
thickness = 240                       # Diamond-layer thickness along the z-direction in nm. The central spin is positioned at the center of the layer.
usp_flag = 0                          # Set usp_flag=1 to simulate P1-centers, usp_flag=0 to simulate electron spins in the bath.
gamma_b = 28000                       # Gyromagnetic ratio of the bath spins in MHz/T. Default is for electron spins.
r_dipole = 45                         # Dipole-radius of the bathspins defining the distance in which dipolar interactions between bath spins are considered.
hf_for_P1 = 1                         # Flag to change the spin bath from electron spins to P1 centres (use 0 for false, 1 for true).
```
Ensure input.py is in the same directory as the main script.

2. **Command-Line Arguments**: You can directly pass simulation parameters when running the script. This method is suitable for quickly changing parameters without modifying the code. For example:

```bash
python main.py --no_systems 10 --number_atoms 100 --t_max 2000 --number_interlaced 5 --max_size 20 --min_size 1 --concentration 1e-6 --thickness 5 --r_dipole 3 --time_step 50 --gamma_b 4 --hf_for_P1 1
```
The available command-line arguments are:
--info, --max_size, --min_size, --no_systems, --t_max, --number_atoms, --number_interlaced, --concentration, --thickness, --r_dipole, --time_step, --gamma_b, --hf_for_P1

## Collaboration
This project is a collaborative effort and contributions are welcome. To participate or inquire, please contact the project team.

## Acknowledgments
We thank Mrunal Divecha for his help with restructuring and optimizing the code. 

## License
This project is licensed under the BSD 3-Clause License. For more details, see the LICENSE file in the project repository.



