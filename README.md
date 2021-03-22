# get-binding-pockets-dataset
A repo with the source code to get a dataset of real binding sites composed by atoms distances matrices and atoms types vectors.

To create the dataset you should first run `python get_nearest_atoms_from_ligand.py` that will extract the nearer protein atoms to ligands in the complexes present in the *PDB_BindingDB.tsv* database. Then  you shoul run `python get_binding_site_distance_matrix_types_vector.py` to get the actual dataset.
