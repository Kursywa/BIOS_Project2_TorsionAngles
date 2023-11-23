# BIOS_Project2_TorsionAngles
## Required packages
The script should be launched in the Linux terminal. If using the Windows operating system, it is advised to use WSL with these required packages installed:
* python 3.1.6
* Biopython 
* matplotlib
* dssp
* math
* os
* sys
## Setup
1. Clone the repository on your local machine
2. In linux terminal navigate to project directory
3. Run the script via following command:
   ```bash
   python3 Torsion_angles.py [mode] [pdb_id]
   ```
note that mode can be either "protein" or "rna" depending on what type of structure is going to be analyzed. 
pdb_id can be any 4 characters long string that represents a RNA or protein structure identifier
4. If there is need to analyze other protein or RNA structure download it from https://www.rcsb.org/ and put the file in working directory

