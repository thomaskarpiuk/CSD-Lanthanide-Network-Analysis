# About
Three scripts for analyzing the geometries of lanthanide complexes in the Cambridge Structural Database. Initially built for analyzing eight-coordinate centres (LnX8).

1. CSD_LnX8_extractor. Extracts eight-coordinate lanthanide centres from the Cambridge Structural Database (CSD), filters out incorrect, disordered, or duplicate hits, and writes the atom coordinates to file (.cor) for analysis with Shape v2.1 software.
2. mol_fragment_extractor. Takes structures from the CSD and creates a .mol2 file containing fragments of all atoms within a certain number of bonds (shells) from that central atom.
3. similarity_matrix_constructor. Constructs a similarity matrix-based network of nodes and edges, based on similarity scores between structure fragments.

These programs are associated with a paper on geometry trends of lanthanide complexes in the Cambridge Structural Database. Details will follow upon publication.

For inquiries, please email at tkarpiuk@sfu.ca
