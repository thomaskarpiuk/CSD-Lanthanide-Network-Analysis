#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Takes a .mol2 file containing molecular fragments and outputs data for a node network based on similarity between
the fragments; outputs the node and edge information based on a chosen similarity threshold.

@author: Thomas Karpiuk
Department of Chemistry, Simon Fraser University
e-mail: tkarpiuk@sfu.ca
"""

import time
import logging
from ccdc.io import MoleculeReader
from ccdc.search import SimilaritySearch

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------
SIM_THRESHOLD = 0.75
BUFFER_WINDOW = 5
MIN_NUM_ATOMS = 9
INPUT_FILE_NAME = "LnX8_frag_molecules.mol2"
OUTPUT_FILE_NAME = "LnX8_similarity_matrix"


# ------------------------------------------------------------------------
# Function 'sort molecules by number of atoms'
# ------------------------------------------------------------------------
def sort_molecules_num_atoms(input_file):
    """
    Takes a .mol2 file and sorts the molecule objects into a dictionary based on
    the number of atoms in the molecule

    Parameters:
        input_file (str): the name of the .mol2 input file.
    
    Returns:
        atom_count_dict: dictionary containing molecule objects sorted by number of atoms.
        max_num_atoms: integer of the largest number of atoms encountered
    """

    atom_count_dict = {}
    max_num_atoms = 0
    for mol in MoleculeReader(input_file):

        # Add mol to the atom count dictionary based on size (number of atoms)
        num_atoms = len(mol.atoms)
        if num_atoms in atom_count_dict:
            atom_count_dict[num_atoms].append(mol)
        else:
            atom_count_dict[num_atoms] = [mol]

        # Update max number of atoms
        if num_atoms > max_num_atoms: 
            max_num_atoms = num_atoms
    
    return atom_count_dict, max_num_atoms


# ------------------------------------------------------------------------
# Function 'merge duplicates'
# ------------------------------------------------------------------------
def merge_duplicates(output_file, atom_count_dict):
    """
    Takes the atom_count_dict and merges all duplicate molecules, which gives the nodes of the network. ALso writes the
    duplicate info to file.

    Parameters:
        output_file (str): the name of the output file for logging the duplicate info
        atom_count_dict (dict): dictionary containing molecule objects sorted by number of atoms.
    
    Returns:
        reduced_atom_count_dict: dictionary containing molecule objects sorted by number of atoms,
        with duplicate molecules merged
        node_dict: dictionary containing all of the nodes: a tuple of identifier, Ln label and the integer size of the node
    """

    reduced_atom_count_dict = {}
    node_dict = {}
    # Search within each molecule size for duplicates
    with open(f"{output_file}_group_info.csv", "w") as f:
        for key in atom_count_dict.keys():
            n = len(atom_count_dict[key])
            search_mol_list = atom_count_dict[key]
            for i in range(n):
                if not search_mol_list: continue

                # Remove the search_mol from the list and add it to the reduced dict
                search_mol = search_mol_list.pop(0)
                if key in reduced_atom_count_dict:
                        reduced_atom_count_dict[key].append(search_mol)
                else:
                    reduced_atom_count_dict[key] = [search_mol]
                
                # If search_mol is the last mol left, add it to the node_dict as a unique entry
                if not search_mol_list:
                    node_dict[(search_mol.identifier, search_mol.atoms[0].label)] = 1
                    continue

                # Start the search
                duplicate_query = SimilaritySearch(search_mol)
                duplicate_query.threshold = 1.0
                duplicate_hits = duplicate_query.search(database=search_mol_list, max_hit_structures=None, max_hits_per_structure=None)

                # If there are no duplicates, add the search_mol to the node_dict as a unique entry
                if not duplicate_hits:
                    node_dict[(search_mol.identifier, search_mol.atoms[0].label)] = 1
                    continue

                # If there are duplicates, print duplicate info to file and then add just the parent molecule (search_mol) to the node_dict with the group size
                group_size = len(duplicate_hits) + 1
                group_list = ["_".join((search_mol.identifier, search_mol.atoms[0].label))]
                for duplicate in duplicate_hits:
                    group_list.append("_".join((duplicate.identifier, duplicate.molecule.atoms[0].label)))

                    # Remove the duplicate molecules from the search_mol_list
                    search_mol_list.remove(duplicate.molecule)
                    
                print(f"{group_size},{','.join(group_list)}", file=f)
                node_dict[(search_mol.identifier, search_mol.atoms[0].label)] = group_size
    
    return reduced_atom_count_dict, node_dict


# ------------------------------------------------------------------------
# Function 'construct_sim_matrix'
# ------------------------------------------------------------------------
def construct_sim_matrix(output_file, sim_threshold, buffer_window, min_num_atoms, reduced_atom_count_dict, max_num_atoms):
    """
    Takes the reduced_atom_count_dict and calculates the similarity between all molecules (nodes) within a range of number of atoms.
    Writes all similarities between nodes above a threshold to file as the edges of the network.

    Parameters:
        output_file (str): the name of the output file for logging the duplicate info
        sim_threshold (float): the threshold similarity score to use for the network
        buffer_window (int): the number of (+/-) atoms to search against for molecules
        min_num_atoms (int): the minimum number of atoms expected for molecules in the dataset, e.g., 9 in the case of LnX8 fragments.
        reduced_atom_count_dict (dict): dictionary containing unique molecule objects (representing the nodes of the network)
        sorted by number of atoms.
        max_num_atoms (int): the maximum number of atoms calculated for molecules in the dataset
    
    Returns:
        reduced_atom_count_dict: dictionary containing molecule objects sorted by number of atoms,
        with duplicate molecules merged
        node_dict: dictionary containing all of the nodes: a tuple of identifier, Ln label and the integer size of the node
    """

    buffered_node_dict = {}
    for num_atoms in range(min_num_atoms, max_num_atoms-buffer_window):
        if num_atoms not in reduced_atom_count_dict: continue

        # Base group is the molecules with the size that you are searching
        base_group_len = len(reduced_atom_count_dict[num_atoms])
        key = (num_atoms, base_group_len)
        buffered_node_dict[key] = []

        # Add in the molecules with larger size than the base group, according to the buffer window
        for i in range(num_atoms, num_atoms + buffer_window):
            if i not in reduced_atom_count_dict: continue
            buffered_node_dict[key].extend(reduced_atom_count_dict[i])

    # Prepare for the similarity search
    with open(f"{output_file}_edges.csv", "w") as f:

        # Search each base group by similarity and write to file
        for key in buffered_node_dict.keys():
            (num_atoms, base_group_len) = key
            search_mol_list = buffered_node_dict[key]

            # Searching only the molecules in the base group
            for i in range(base_group_len):
                if not search_mol_list: continue
                search_mol = search_mol_list.pop(0)
                if not search_mol_list: continue

                # Start the search
                sim_query = SimilaritySearch(search_mol)
                sim_query.threshold = sim_threshold
                sim_hits = sim_query.search(database=search_mol_list, max_hit_structures=None, max_hits_per_structure=None)

                if not sim_hits: continue

                # If not unique, add each hit to the output
                for h in sim_hits:
                    print(f"{search_mol.identifier}_{search_mol.atoms[0].label},{h.identifier}_{h.molecule.atoms[0].label},{h.similarity}", file=f)


def main():
    start_time = time.time()

    # ------------------------------------------------------------------------
    # Read the molecules from file, merge duplicates, and sort them into a dictionary by number of atoms
    # ------------------------------------------------------------------------
    logging.info("Sorting molecules...")
    atom_count_dict, max_num_atoms = sort_molecules_num_atoms(input_file=INPUT_FILE_NAME)

    logging.info(f"Merging duplicates...")
    reduced_atom_count_dict, node_dict = merge_duplicates(output_file=OUTPUT_FILE_NAME, atom_count_dict=atom_count_dict)
    
    # ------------------------------------------------------------------------
    # Print the node info to file
    # ------------------------------------------------------------------------
    logging.info(f"Print node info to file {OUTPUT_FILE_NAME}_all_nodes.csv...")
    with open(f"{OUTPUT_FILE_NAME}_all_nodes.csv", "w") as f:
        for (identifier, label), group_size in node_dict.items():
            print(f"{identifier}_{label},{group_size}", file=f)

    # ------------------------------------------------------------------------
    # Construct the similarity matrix network and print the edges info to file
    # ------------------------------------------------------------------------
    logging.info("Running the similarity search...")
    construct_sim_matrix(OUTPUT_FILE_NAME, SIM_THRESHOLD, BUFFER_WINDOW, MIN_NUM_ATOMS, reduced_atom_count_dict, max_num_atoms)

    elapsed_time = time.time() - start_time
    logging.info(f"Execution time: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    main()
