#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Given a list of tuples of identifier and central atom, create a mol2 file containing fragments of all atoms
within a certain number of bonds (shells) from that central atom.

@author: Thomas Karpiuk
Department of Chemistry, Simon Fraser University
e-mail: tkarpiuk@sfu.ca
"""

import json
import time
import logging
from ccdc.io import MoleculeWriter, MoleculeReader
from ccdc.molecule import Molecule

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------
STRUCTURE_LIST_FILENAME = "LnX8_hitlist.json"
STRUCTURE_LIST_NAME = "LnX8_hitlist"
OUTPUT_FILE_NAME = "LnX8_frag_molecules.mol2"
NUMBER_OF_SHELLS = 3


# ------------------------------------------------------------------------
# Function 'open structure list'
# ------------------------------------------------------------------------
def open_structure_list(structure_list_filename, structure_list_name):
    """
    Load a .json file containing a list of structures.

    Parameters:
        structure_list_filename (str): the .json file that contains the list of structures.
        structure_list_name (str): the name of the structure list, which is a list of tuples of identifier and central atom label for each hit
    
    Returns:
        structure_list: a list of tuples of identifier and central atom label for each structure
    """
    with open(structure_list_filename, "r") as f:
        data = json.load(f)
    structure_list = [tuple(pair) for pair in data[structure_list_name]]

    return structure_list


# ------------------------------------------------------------------------
# Function 'add shells'
# ------------------------------------------------------------------------
def add_shells(frag_mol, mol, central_atom, max_shells=5):
    """
    Add atoms and bonds to a fragment molecule up to a given shell depth.

    Parameters:
        frag_mol (Molecule): The fragment molecule to which atoms and bonds will be added.
        mol (Molecule): The CSD molecule to make reference to.
        central_atom (Atom): The central atom from which the shells will be built.
        max_shells (int): The maximum number of shells to add.
    """

    # Add the central atom to the fragment molecule first
    frag_central_atom = frag_mol.add_atom(central_atom)
    
    # Use a set to keep track of atoms added to the fragment molecule
    frag_mol_set = set()
    frag_mol_set.add(central_atom)

    # Keep track of atoms added to each shell, as them in the CSD entry and in the newly built fragment
    shell_atoms = [[(central_atom, frag_central_atom)]]

    # Iterate over each shell, up to the max_shells depth
    for shell in range(1, max_shells + 1):
        new_shell_atoms = []

        # For each atom in the previous shell
        for (base_atom, frag_base_atom) in shell_atoms[shell - 1]:
            next_shell = base_atom.neighbours
            added_atoms = []

            # For each atom in the next shell
            for next_atom in next_shell:

                # Check if the atom has already been added to the fragment
                if next_atom not in frag_mol_set:

                    # if not, add the atom to the fragment
                    frag_next_atom = frag_mol.add_atom(next_atom)
                    frag_mol_set.add(next_atom)
                    
                    # Add a bond between the atom and its parent atom in the previous shell
                    try:
                        mol_bond = mol.bond(base_atom.label, next_atom.label).bond_type
                    except RuntimeError:
                        mol_bond = 'Single'
                    frag_mol.add_bond(mol_bond, frag_base_atom, frag_next_atom)

                    added_atoms.append((next_atom, frag_next_atom))
                else:
                    # If the atom is already in the fragment, check to add any missing bonds
                    try:
                        frag_next_atom = frag_mol.atom(next_atom.label)
                        frag_mol.add_bond('Single', frag_base_atom, frag_next_atom)
                    except RuntimeError:
                        continue

            # Add the new atoms to the next shell
            new_shell_atoms.extend(added_atoms)

        # Stop if no new atoms were added in the current shell
        if not new_shell_atoms:
            break

        shell_atoms.append(new_shell_atoms)

def main():
    start_time = time.time()
    counter = 0
    mol_reader = MoleculeReader('CSD')

    # ------------------------------------------------------------------------
    # Load a structure list
    # ------------------------------------------------------------------------
    logging.info("Opening structure list...")
    structure_list = open_structure_list(STRUCTURE_LIST_FILENAME, STRUCTURE_LIST_NAME)

    # ------------------------------------------------------------------------
    # Build a fragment for each structure
    # ------------------------------------------------------------------------
    with MoleculeWriter(OUTPUT_FILE_NAME, append=True) as writer:
        for structure in structure_list:
            (i, Ln) = structure

            # Define the crystal, molecule, and Ln atom
            mol = mol_reader.molecule(i)
            mol.remove_hydrogens()
            try:
                Ln_atom = mol.atom(Ln)
            except RuntimeError:
                for atom in mol.atoms:
                    if atom.label == Ln and len(atom.neighbours) == 8:
                        Ln_atom = atom
            frag_mol = Molecule(identifier=(i))

            # Build the fragment molecule and add to the fragment molecule list
            add_shells(frag_mol, mol, Ln_atom, max_shells=NUMBER_OF_SHELLS) 
            frag_mol.atoms[0].atomic_symbol = 'Ac' # Using Ac as a dummy atom to replace the central atom
            writer.write(frag_mol)

            counter += 1
            if counter % 1000 == 0:
                logging.info(f"Processed {counter:,} structures...")

    elapsed_time = time.time() - start_time
    logging.info(f"Execution time: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    main()