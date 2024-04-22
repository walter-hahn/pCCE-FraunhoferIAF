
# Partition-based Cluster-Correlation Expansion pCCE

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
1. **Command-Line Arguments**: You can directly pass simulation parameters when running the script. This method is suitable for quickly changing parameters without modifying the code. For example:

```bash
python main.py --no_systems 10 --number_atoms 100 --t_max 2000 --number_interlaced 5 --max_size 20 --min_size 1 --base_concentrations [1e-6, 2e-6] --thickness 5 --r_dipole 3 --time_step 50 --gamma_b 4 --hf_for_P1 1
```
The available command-line arguments include:

--info: Display information about input methods.

--max_size: Maximum cluster size.

--min_size: Minimum cluster size.

--no_systems: Number of systems.

--t_max: Maximum time for simulation.

--number_atoms: Number of atoms in a system.

--number_interlaced: Number of interlaced simulations.

--base_concentrations: Base concentrations, specify multiple values separated by spaces.

--thickness: Layer thickness of the spin distribution.

--r_dipole: Dipole weight.

--time_step: Time step in milliseconds.

--gamma_b: Gamma B value.

--hf_for_P1: Flag to change the spin bath from Electron spins to P1 centres (use 0 for false, 1 for true).

2. **Modifying input.py**: For more persistent changes, you can modify the input.py file with your desired parameters. This file overrides the default settings specified in constants.py. If input.py does not exist, the simulation will fall back to using the default values.

Sample input.py configuration:
```bash
max_part_size = 4                     # Maximum partition size    
min_part_size = 4                     # Minimum partition size

number_systems = 1                    # Number of systems considered. Each system corresponds to a random spatial distribution of spins. MPI parallelization distributes theses systems among available cores.
t_max = 500                           # Maximum time in µs
number_internal_avg = 20              # Number of averages performed internally (Monte Carlo bath state sampling)
number_spins = 144                    # Total number of bath-spins included in the calculation       
base_concentrations = [1e-6, 2e-6]    # List of concentrations (Fraction of atoms replaced by defects)
thickness = 240                       # Diamond-layer thickness along the z-direction in nm. The central spin is positioned at the center of the layer.
usp_flag = 0                          # Set usp_flag=1 to simulate P1-centers, usp_flag=0 to simulate electron spins in the bath.
gamma_b = 28000                       # Gyromagnetic ratio of the bath spins in MHz/T. Default is for electron spins.
r_dipole = 45                         # Dipole-radius of the bathspins defining the distance in which dipolar interactions between bath spins are considered.
```
Ensure input.py is in the same directory as the main script.

## Collaboration
This project is a collaborative effort, and contributions are welcome. To participate or inquire, please contact the project team.

## Acknowledgments
We thank Mrunal Divecha for helping with restructuring and optimizing the code. 

## Caution
The code is under continuous development and, hence, comes without guarantee.

## License
This project is licensed under the BSD 3-Clause License. For more details, see the LICENSE file in the project repository.



