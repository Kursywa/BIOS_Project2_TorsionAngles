import sys
from Bio import PDB
import numpy as np
from Bio.PDB import PDBList
from Bio.PDB.vectors import calc_dihedral
from math import degrees
from Bio.PDB.DSSP import DSSP
import matplotlib.pyplot as plt
from Bio.PDB import Polypeptide
import os

def main():
    # checking arguments in terminal
    mode, pdb_id = check_command_input()

    # Downloading the wanted structure form PDB database
        # obsolete cause pdblist gave weird format aso dssp couldnt read it
    # pdb_filename = download_structure(pdb_id)
    pdb_filename = "./" + pdb_id + ".pdb"

    print(pdb_filename)
    # Parse the structure from the downloaded file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_filename)

    if mode == "rna":
        # Calculate torsion angles
        residue_angles = process_structure(structure)

        # Print the results to output
        print_test(residue_angles)

        # Save data to matrix
        matrix = save_matrix(residue_angles, structure)

        # Printing the matrix for testing
        print(matrix)

        # Saving results for angles in a file
        save_to_file(matrix,"RNA_angles.txt")

    elif mode == "protein": # mode configured at the start
        phi_angles, psi_angles, secondary_structure = calculate_phi_psi(structure, pdb_filename)

        filename, file_extension = os.path.splitext(os.path.basename(pdb_filename))
        plot_ramachandran(phi_angles, psi_angles, secondary_structure,
                          os.path.basename(filename + "_" + file_extension[1:]))


def check_command_input():# get command line arguments 1- mode either protein/RNA 2- pdb ID
    if len(sys.argv) != 3:
        print("Usage: python Torsion_Angles.py protein/RNA [pdb_ID]")
        sys.exit(1)
    elif sys.argv[1] not in ["protein", "rna"]:
        print("Error: Argument must be 'protein' or 'RNA'.")
        sys.exit(1)
    elif len(sys.argv[2]) != 4:
        print("Error: PDB identifier should be 4 characters long")
        sys.exit(1)
    else:
        mode = sys.argv[1]
        pdb_id = sys.argv[2]
    return mode, pdb_id

def download_structure(pdb_id):
    # Download the structure using the PDB ID
    pdbl = PDBList() # necessary structure
    pdb_filename = pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=".") # Obtaining specific pdb file
    return pdb_filename


def process_structure(structure):
    angles = { # exact atoms which are taken into account in calculating torsion angles
        "alpha": ("O3'", "P", "O5'", "C5'"),
        "beta": ("P", "O5'", "C5'", "C4'"),
        "gamma": ("O5'", "C5'", "C4'", "C3'"),
        "delta": ("C5'", "C4'", "C3'", "O3'"),
        "epsilon": ("C4'", "C3'", "O3'", "P"),
        "zeta": ("C3'", "O3'", "P", "O5'"),
        "chi_purines": ("O4'", "C1'", "N9", "C4"),
        "chi_pyrimidines": ("O4'", "C1'", "N1", "C2"),
    }

    all_residue_angles = {}

    for model in structure:
        for chain in model:
            for i, residue in enumerate(chain):
                angles_values = {}
                #From the dictionary above
                for angle_name, atoms in angles.items():
                    try:
                        match angle_name:
                            # cases are for every angle that has to take into account previous or next residue
                            case "alpha":
                                if (i == 0):
                                    raise NoPrev
                                elif(i > 0):
                                    previous_residue = chain.get_list()[i - 1]
                                    atoms_objects = [previous_residue[atoms[0]]] + [residue[atoms[i]] for i in range(1, 4)]
                            case "epsilon":
                                next_residue = chain.get_list()[i + 1]
                                atoms_objects =  [residue[atoms[i]] for i in range(0, 3)] + [next_residue[atoms[3]]]
                            case "zeta":
                                next_residue = chain.get_list()[i + 1]
                                atoms_objects = [residue[atoms[i]] for i in range(0, 2)] + [next_residue[atoms[i]] for i in range(2,4)]
                            case _: # every other angle: beta, gamma, delta, chi
                                atoms_objects = [residue[atoms[i]] for i in range(4)]
                        # obtaining vectors of said atoms
                        vectors = [atom.get_vector() for atom in atoms_objects]
                        angle = calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3])
                        # Calc_dihedral returns radians - converting it to degrees
                        angle_degrees = round(degrees(angle), 1)
                        # Updating the temp dictionary
                        angles_values[angle_name] = angle_degrees
                    except Exception as NoPrev:
                        angles_values[angle_name] = "-"
                # Saving only the residues that have at least one angle
                if any(angle != "-" for angle in angles_values.values()):
                    all_residue_angles[residue.id[1]] = angles_values
    return all_residue_angles

def print_test(residue_angles):
    print("res no.\talpha\tbeta\tgamma\tdelta\tepsilon\tzeta\tχ_purines\tχ_pyrimidines")
    for residue_num, angles_values in residue_angles.items():
        print(f"{residue_num}\t"
            f"{angles_values.get('alpha')}\t" # print either angle or null in form of '-'
            f"{angles_values.get('beta')}\t"
            f"{angles_values.get('gamma')}\t"
            f"{angles_values.get('delta')}\t"
            f"{angles_values.get('epsilon')}\t"
            f"{angles_values.get('zeta')}\t"
            f"{angles_values.get('chi_purines')}\t"
            f"{angles_values.get('chi_pyrimidines')}")

def save_matrix(all_residue_angles, structure):
    matrix = np.empty((len(all_residue_angles), 7), dtype = object) # 7 - hardcoded number of examined angles
    # Dtype cause we dont want conversion errors from string to float
    for i, residue in enumerate(structure.get_residues()):
        index = i
        # Not taking more residues than actually calculated
        if index == len(all_residue_angles):
            break
        angles_values = all_residue_angles.get(residue.id[1], {})
        matrix[index, 0] = angles_values.get("alpha", "-")
        matrix[index, 1] = angles_values.get("beta", "-")
        matrix[index, 2] = angles_values.get("gamma", "-")
        matrix[index, 3] = angles_values.get("delta", "-")
        matrix[index, 4] = angles_values.get("epsilon", "-")
        matrix[index, 5] = angles_values.get("zeta", "-")

        # exception for chi variation - displaying only the the correctly associated one
        if any(letter in residue.get_resname() for letter in ["A", "G"]):
            matrix[index, 6] = angles_values.get("chi_purines", "-")
        elif any(letter in residue.get_resname() for letter in ["C", "U"]):
            matrix[index, 6] = angles_values.get("chi_pyrimidines", "-")
        else:
            matrix[index, 6] = "-"
    return matrix

def save_to_file(matrix, output_file_path):
    np.savetxt(output_file_path, matrix, fmt='%s', delimiter='\t', header="α\tβ\tγ\tδ\tε\tζ\tχ", comments='',
               encoding='utf-8')

def calculate_phi_psi(structure, pdb_filename):
    # List initilization
    phi_angles = []
    psi_angles = []
    secondary_structure = []

    for model in structure:
        for chain in model:
            polypeptides = Polypeptide.Polypeptide(
                chain)  # creating Polupeptide object for chain for it to be analyzed
            phi_psi_angles = polypeptides.get_phi_psi_list()  # getting phi psi angles
            dssp = PDB.DSSP(model, pdb_filename)

            residues = list(chain.get_residues())  # getting residues
            for phi, psi in phi_psi_angles:
                if phi is not None and psi is not None:
                    residue_index = polypeptides.get_phi_psi_list().index(
                        (phi, psi))  # getting residue index via torsion angles
                    residue_id = (chain.id, residues[residue_index].get_id()[
                        1])  # getting residue identifier made od chain_id and residue index
                    phi_angles.append(np.degrees(phi))  # putting phi and psi angles into list converting them
                    psi_angles.append(np.degrees(psi)) # a priori to degrees from radians
                    secondary_structure.append(
                        dssp[residue_id][2])  # dssp used for secondary structure

    return phi_angles, psi_angles, secondary_structure

def plot_ramachandran(phi_angles, psi_angles, secondary_structure, pdb_filename):
    plt.figure(figsize=(8, 6))  # creating new figure object 8x6 inches

    # Mapping all types of secondary structures but highlighting helices as red and beta sheets as blue
    color_map = {'H': 'red', 'E': 'blue'}

    # Rscattering points on the plot and coloring them according to colormap
    for phi, psi, ss in zip(phi_angles, psi_angles, secondary_structure):
        plt.scatter(phi, psi, s=10, alpha=0.5, edgecolors='none', c=color_map[ss])  # Adding points to plot

    # Creating plot annotations
    plt.title('Ramachandran Plot with Secondary Structure Coloring')
    plt.xlabel('Phi Angle (degrees)')
    plt.ylabel('Psi Angle (degrees)')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    legend_labels = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label)
                     # Creating
                     for label, color in color_map.items()]
    plt.legend(handles=legend_labels, title='Secondary Structure',
               loc='upper right')  # Creating legend
    plt.savefig(pdb_filename + '_ramachandran_plot.pdf')  # Saving to file. File should not be in use
    plt.show()

if __name__ == "__main__":
    main()
