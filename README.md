
# Quantum Spin Simulation at Fraunhofer IAF

## Project Overview

This project at Fraunhofer IAF, delves into the complex realm of quantum spin dynamics. The primary focus is on the Hahn echo problem, a phenomenon significant in the study of quantum coherence and decoherence in spin systems. Utilizing advanced computational methods, the project simulates the interactions within spin clusters, exploring how these micro-level interactions influence the macroscopic properties of the system. Supervised by Walter Hahn and in collaboration with Philip Schaetzle from the University of Freiburg, the project stands at the forefront of quantum computational research, aiming to deepen our understanding of quantum mechanics and its potential applications in technology.
## Features

- Creation of Hamiltonian matrices for spin clusters.
- Time propagation of quantum states.
- Calculation of probabilities for Hahn echo experiments.
- Customizable simulation parameters via `constants.py`.

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
python main.py --no_systems 10 --number_atoms 100 --t_max 2000 --number_interlaced 5 --max_size 20
```
Replace the values as per your requirements. The available arguments are --no_systems, --number_atoms, --t_max, --number_interlaced, and --max_size.

2. **Modifying input.py**: For more persistent changes, you can modify the input.py file with your desired parameters. This file overrides the default settings specified in constants.py. If input.py does not exist, the simulation will fall back to using the default values.

Sample input.py configuration:
```bash
no_systems = 10
number_atoms = 100
t_max = 2000
number_interlaced = 5
max_size = 20
```
Ensure input.py is in the same directory as the main script.

## Collaboration
This project is a collaborative effort, and contributions are welcome. To participate or inquire, please contact the project team.

## Acknowledgments
Special thanks to Walter Hahn and Philip Schaetzle for their invaluable guidance and collaboration in this project.

## License
This project is licensed under the GNU Lesser General Public License v2.1. For more details, see the LICENSE file in the project repository.



