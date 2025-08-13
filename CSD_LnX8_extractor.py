#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extracts eight-coordinate lanthanide centres from the Cambridge Structural Database (CSD), filters
out incorrect, disordered, or duplicate hits, and writes the atom coordinates to file for analysis
with Shape v2.1 software.

@author: Thomas Karpiuk
Department of Chemistry, Simon Fraser University
e-mail: tkarpiuk@sfu.ca
"""

import json
import time
import logging
from ccdc.search import SubstructureSearch, QuerySubstructure, QueryAtom, QueryBond

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------
MAX_R_FACTOR = 10
LN_UEQ_THRESHOLD = 0.070
OUTPUT_FILENAME = "LnX8_hitlist.json"
SHAPE_FILENAME = "LnX8_hitlist.cor"

# Lists of atomic symbols for search criteria
LN_ATOMS = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
LIGAND_ATOMS = ['N', 'O', 'F', 'P', 'S', 'Cl', 'Br']
CH_ATOMS = ['B', 'C', 'N', 'O', 'Si', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Te', 'I']

# List of specific problem structures to exclude from analysis
EXCLUDE_LIST = {
    ('Pr1', 'IPICAN'), ('Nd1', 'IPICOB'), ('Sm1', 'IPIDES'),
    ('Eu1', 'IPIDIW'), ('Tb1', 'IPIFAQ'), ('Er1', 'IPIFEU'),
    ('Dy1', 'IPIFOE'), ('Lu1', 'IPIGAR'), ('Yb1', 'IPIGOF'),
    ('Ho1', 'IPIGUL'), ('Tm1', 'IPIHEW'), ('Tb1', 'JUHFOJ'),
    ('Dy1', 'JUHFUP'), ('Ho1', 'JUHGAW'), ('Er1', 'JUHGEA'),
    ('Tm1', 'JUHGIE'), ('Yb1', 'JUHGOK'), ('Dy1', 'NOBRII'),
    ('Tb1', 'XUCYIF'), ('Dy1', 'XUCYOL'), ('Ho1', 'XUCYUR'),
    ('Tm1', 'XUCZEC'), ('Gd2', 'ZAQBUS'), ('Yb1', 'FAKHAB'),
    ('Nd1', 'ATEBIM'), ('Ho3', 'BIQNUN'), ('Er3', 'BIQPAV'),
    ('Sm2', 'XAYNUJ'), ('Eu1', 'INANUJ')
}

# Chemical names that indicate an uncommon oxidation state (non-III) of the lanthanide
INVALID_CHEMICAL_NAMES = {'cerium(iv)', 'terbium(iv)', 'europium(ii)', 'ytterbium(ii)'}

# ------------------------------------------------------------------------
# Function 'create substructure search'
# ------------------------------------------------------------------------
def create_substructure_search(central_atom_list, ligand_atom_list, num_bonds=8, chelate_atom_list=None):
    """
    Creates a substructure search for n-coordinate metal centres, optionally including
    one chelate atom to search for hits with 4-membered chelate rings, which is useful
    to filter out such hits from the parent dataset. Hits with 3-membered chelate rings
    are sorted out be default.

    Parameters:
        central_atom (list): List of atomic symbols for the central atom.
        ligand_atoms (list): List of atomic symbols for the ligand atoms.
        num_bonds (int): coordination number of the central atom.
        chelate_atoms (list, optional): atomic symbols for the chelate ring atom.
    
    Returns:
        SubstructureSearch: Configured substructure search object.
    """
    search = SubstructureSearch()
    substructure = QuerySubstructure()

    # Define substructure atoms
    central_atom = QueryAtom(central_atom_list)
    central_atom.num_bonds = num_bonds
    ligands = [QueryAtom(ligand_atom_list) for _ in range(num_bonds)]
    
    # Add atoms to the substructure
    substructure.add_atom(central_atom)
    for ligand in ligands:
        substructure.add_atom(ligand)

    # Add bonds to the substructure
    for ligand in ligands:
        bond = QueryBond('Single')
        bond.bond_smallest_ring = ('!=', 3)
        substructure.add_bond(bond, central_atom, ligand)

    # If adding the four-membered chelate ring, add the chelate atom and bonds to the substructure
    if chelate_atom_list:
        chelate = QueryAtom(chelate_atom_list)
        substructure.add_atom(chelate)
        for i in range(2):
            substructure.add_bond(QueryBond(), ligands[i], chelate)

    # Add the substructure to the search and configure settings
    search.add_substructure(substructure)
    search.settings.has_3d_coordinates = True
    search.settings.max_r_factor = MAX_R_FACTOR
    return search

# ------------------------------------------------------------------------
# Function 'filter hits'
# ------------------------------------------------------------------------
def filter_hits(hits, chelate_ring_hits, filename):
    """
    Filters hits based on multiple criteria and writes each successful hit to .cor file formatted
    for input into the Shape v2.1 software.

    Parameters:
        hits (list): List of :class:`ccdc.search.SubstructureSearch.SubstructureHit` instances.
        chelate_ring_hits (list): List of tuples containing identifier and Ln atom label for 4-ch search hits.
        filename (str): Name of the output .cor file.
    
    Returns:
        filtered_hits: a list of tuples of identifier and Ln atom label for filtered hits.
    """
    counter = 0
    filtered_hits = []
    processed_ids = {}
    base_id_hits = {}
    statistics = {
        "exclude_list": 0, "chelate_rings": 0, "occupancy": 0,
        "ueq": 0, "duplicate_sg": 0, "chemical_name": 0
    }
    
    with open(filename, "w") as f:
        for h in hits:
            counter += 1
            if (counter % 1000 == 0):
                logging.info(f"at hit {counter}, {len(filtered_hits)} added to filtered hits so far...")
            Ln_atom = h.match_atoms()[0]
            
            # Remove specific problem structures
            if (Ln_atom.label, h.identifier) in EXCLUDE_LIST:
                statistics["exclude_list"] += 1
                continue
            
            # Remove hits that contain 4-membered chelate rings
            if (h.identifier, Ln_atom.label) in chelate_ring_hits:
                statistics["chelate_rings"] += 1
                continue
            
            # Remove hits where any of the atoms do not have occupancy of 1 (disordered)
            if any(a.occupancy != 1 for a in h.match_atoms()):
                statistics["occupancy"] += 1
                continue
            
            # Remove hits with too large displacement of the Ln atom (likely disordered)
            Ln_ueq = getattr(Ln_atom.displacement_parameters, "isotropic_equivalent", None)
            if Ln_ueq and Ln_ueq >= LN_UEQ_THRESHOLD:
                statistics["ueq"] += 1
                continue
            
            # Remove hits that are duplicates; i.e, with the same identifier and space group
            base_id = h.identifier[:6]
            current_sg = h.crystal.spacegroup_symbol
            if base_id in processed_ids:
                base_sg = processed_ids[base_id]
                if len(h.identifier) > 6 and current_sg == base_sg:
                    statistics["duplicate_sg"] += 1
                    continue
            
            # Remove hits that are not in the common (iii) oxidation state
            chemical_name = h.entry.chemical_name
            if chemical_name and any(kw in chemical_name.lower() for kw in INVALID_CHEMICAL_NAMES):
                statistics["chemical_name"] += 1
                continue
            
            # If its a new structure, store its space group and initialize a counter for that identifier
            if base_id not in processed_ids:
                processed_ids[base_id] = current_sg
                base_id_hits[base_id] = 0
            base_id_hits[base_id] += 1

            # Append the identifier, Ln atom label, and hit number for that identifier
            filtered_hits.append((h.identifier, Ln_atom.label))

            # Write the hit details to the .cor file for Shape v2.1 analysis
            print(f"{h.identifier:<10}  **FRAG**        {base_id_hits[base_id]}", file=f)

            # Print the Ln atom label and coordinates
            print(f"{Ln_atom.label:<10}  {Ln_atom.coordinates.x:10.5f}  {Ln_atom.coordinates.y:10.5f}  {Ln_atom.coordinates.z:10.5f}", file=f)
            
            # Print the eight ligand atom labels and coordinates
            for a in Ln_atom.neighbours:
                print(f"{a.label:<10}  {a.coordinates.x:10.5f}  {a.coordinates.y:10.5f}  {a.coordinates.z:10.5f}", file=f)
        
    # Report statistics on the filtering
    for key, value in statistics.items():
        logging.info(f"Removed by {key.replace('_', ' ')}: {value}")
    logging.info(f'number of hits after filtering: {len(filtered_hits)}')
        
    return filtered_hits


def main():
    start_time = time.time()

    # ------------------------------------------------------------------------
    # Perform search for LnX8 hits with any 4-membered chelate rings (to remove later)
    # ------------------------------------------------------------------------
    logging.info("Performing LnX8 4-chelate ring search...")
    LnX8_4ch_search = create_substructure_search(LN_ATOMS, LIGAND_ATOMS, 8, CH_ATOMS)
    hits_4ch = LnX8_4ch_search.search()
    logging.info(f"Hits from LnX8 4-chelate ring search: {len(hits_4ch)}")
    chelate_ring_hits = {(h.identifier, h.match_atoms()[0].label) for h in hits_4ch}
    
    # ------------------------------------------------------------------------
    # Perform LnX8 search
    # ------------------------------------------------------------------------
    logging.info("Performing LnX8 search...")
    LnX8_search = create_substructure_search(LN_ATOMS, LIGAND_ATOMS, 8)
    hits = LnX8_search.search()
    logging.info(f"Initial hits from LnX8 search: {len(hits)}")

    # ------------------------------------------------------------------------
    # Apply filtering and output hits to .cor file for Shape v2.1 software
    # ------------------------------------------------------------------------
    logging.info("Filtering hits...")
    filtered_hits = filter_hits(hits, chelate_ring_hits, SHAPE_FILENAME)

    # ------------------------------------------------------------------------
    # Output dataset as .json file
    # ------------------------------------------------------------------------
    output_data = {
        "LnX8_hitlist": filtered_hits,
        "LnX8_hitlist_ids": list({identifier for identifier, _ in filtered_hits})
    }
    
    try:
        with open(OUTPUT_FILENAME, "w") as f:
            json.dump(output_data, f, indent=4)
        logging.info(f"Data saved to {OUTPUT_FILENAME}")
    except IOError as e:
        logging.error(f"Error saving data: {e}")

    elapsed_time = time.time() - start_time
    logging.info(f"Execution time: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    main()
