# About
Three scripts for analyzing the geometries of lanthanide complexes in the Cambridge Structural Database. Initially built for analyzing eight-coordinate centres (LnX8).

1. CSD_LnX8_extractor. Extracts eight-coordinate lanthanide centres from the Cambridge Structural Database (CSD), filters out incorrect, disordered, or duplicate hits, and writes the atom coordinates to file (.cor) for analysis with Shape v2.1 software.
2. mol_fragment_extractor. Takes structures from the CSD and creates a .mol2 file containing fragments of all atoms within a certain number of bonds (shells) from that central atom.
3. similarity_matrix_constructor. Constructs a similarity matrix-based network of nodes and edges, based on similarity scores between structure fragments.

These programs are associated with the following paper on geometry trends of lanthanide complexes in the Cambridge Structural Database:

Karpiuk, T. E. and Leznoff, D. B. "Strategies to control the geometry and symmetry around lanthanide centres for tailored luminescence and magnetism" Nature Commun. 2026. DOI: https://doi.org/10.1038/s41467-026-69445-6.

These programs were written using the CSD Python API version 3.0.17 and CSD version 5.45 (including up to June 2024 release). The API and database can be installed at the following link: https://www.ccdc.cam.ac.uk/support-and-resources/downloads/

The programs may not work with the latest CSD Python API and CSD database versions -- updates will follow.

For inquiries or suggestions please email me at tkarpiuk@sfu.ca
