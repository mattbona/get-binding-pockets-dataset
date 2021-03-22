import numpy as np
import scipy
from scipy.spatial import distance
from numpy import random
import pandas as pd
import os

def GetRandomNAtomsFromChain(chain_coordinate_matrix, chain_atomtype, natoms):
    selected_atoms = np.random.choice(chain_coordinate_matrix.shape[0], natoms, replace=False)
    random_chain_atoms = chain_coordinate_matrix[selected_atoms, :]
    random_chain_atomtype = chain_atomtype[selected_atoms]

    return random_chain_atoms, random_chain_atomtype

def GetFirstNAtomsFromChain(chain_coordinate_matrix, chain_atomtype, chain_distSqrd, natoms):
    idx_sorted_dists = np.argsort(chain_distSqrd)[:natoms]
    sorted_chain_atoms = chain_coordinate_matrix[idx_sorted_dists, :]
    sorted_chain_atomtype = chain_atomtype[idx_sorted_dists]

    return sorted_chain_atoms, sorted_chain_atomtype


def GetDummyVAtomType(chain_atomtype):
    # recognize types:['C' 'N' 'O' 'S']
    dummy_chain_atomtype = np.array([])
    for atype in chain_atomtype:
        if atype == 'N':
            ndummy = np.array([1, 0, 0, 0, 0])
            dummy_chain_atomtype = np.concatenate((dummy_chain_atomtype, ndummy))
        elif atype == 'C':
            cdummy = np.array([0, 1, 0, 0, 0])
            dummy_chain_atomtype = np.concatenate((dummy_chain_atomtype, cdummy))
        elif atype == 'O':
            odummy = np.array([0, 0, 1, 0, 0])
            dummy_chain_atomtype = np.concatenate((dummy_chain_atomtype, odummy))  
        elif atype == 'S':
            sdummy = np.array([0, 0, 0, 1, 0])
            dummy_chain_atomtype = np.concatenate((dummy_chain_atomtype, sdummy))
        else:
            sdummy = np.array([0, 0, 0, 0, 1])
            dummy_chain_atomtype = np.concatenate((dummy_chain_atomtype, sdummy))
            
    return dummy_chain_atomtype
#################################################################################

dist_size_out = open("distances_matrix_sizes_sort_10.dat","w")
lID_out = open("pdb_chain_ligand_HETID_pname_oname.dat","w")

pdbs_stat = pd.read_csv(
             filepath_or_buffer='pdbs_statistics.dat',
             header=None, 
             sep='____')
             #delim_whitespace=True)
        
print(pdbs_stat)

df = pdbs_stat.drop_duplicates()
pdbs = df.drop_duplicates(subset=0, keep="last")[0].to_numpy()
ligand_HETID = df.drop_duplicates(subset=0, keep="last")[1].to_numpy()
proteinChemName = df.drop_duplicates(subset=0, keep="last")[2].to_numpy()
organismName = df.drop_duplicates(subset=0, keep="last")[3]
#organismName = organismName.str.split(' ').str[0].to_numpy()

unique_protein_name, counts_protein_name = np.unique(proteinChemName, return_counts=True)
unique_org_name, counts_org_name = np.unique(organismName, return_counts=True)

print("Number of unique protein chem names: ", len(unique_protein_name))
print("Number of unique organism names: ", len(unique_org_name))

proteinName_dict = dict(zip(unique_protein_name, counts_protein_name))
orgName_dict = dict(zip(unique_org_name, counts_org_name))
print("Protein chem names: ")
print(np.flip(sorted(proteinName_dict.items(), key = lambda kv: kv[1]))[:10])
print("Organism names: ")
print(np.flip(sorted(orgName_dict.items(), key = lambda kv: kv[1]))[:10])

ligand_HETID = np.vstack((pdbs,ligand_HETID, proteinChemName, organismName)).T

rnd_choice = True
permutation = False
np.random.seed(42)
natoms = 10

pdbs = os.listdir('nearest_chains_pdbs/')
with open('distances_matrix_sort_10.dat','w') as output:
    for ipdb in pdbs:
        print(ipdb)
        ligand_ID = ligand_HETID[np.where(ligand_HETID[:,0] == ipdb[:4]), 1][0][0]
        proteinChemName = ligand_HETID[np.where(ligand_HETID[:,0] == ipdb[:4]), 2][0][0]
        organismName = ligand_HETID[np.where(ligand_HETID[:,0] == ipdb[:4]), 3][0][0]
        print("Computing: ", ipdb[:4], ligand_ID)
        df = pd.read_csv(
             filepath_or_buffer='nearest_chains_pdbs/'+ipdb,
             header=None,
             delim_whitespace=True, #deal with different spaced columns
             skiprows=1, #skip first row
             skipfooter=1) #skip last row

        if df.empty!=True: #if pdb file is not empty...
            #print(df.tail())
            # Print all the chain label
            chains = df[4].unique()
            print("Chains: ", chains)

            # Divide the dataframe by chain label
            dfs = dict(tuple(df.groupby([4])))
            # Iterate over the single chains
            for ichain in chains:
                print(ichain,' chain')

                chain_coordinate_matrix = dfs[ichain].iloc[:,7:10].values
                chain_atomtype = dfs[ichain].iloc[:,2].str[0].values
                chain_distSqrd = dfs[ichain].iloc[:,6].values
                #print(chain_coordinate_matrix.shape, chain_atomtype.shape)
                #print(chain_atomtype)
                if rnd_choice == True: 
                    if chain_coordinate_matrix.shape[0] >= natoms:
                        chain_coordinate_matrix, atomtype_arrays = GetFirstNAtomsFromChain(chain_coordinate_matrix,
                                                                                           chain_atomtype, chain_distSqrd,
                                                                                           natoms)
                        #chain_coordinate_matrix, atomtype_arrays = GetRandomNAtomsFromChain(chain_coordinate_matrix, 
                        #                                                                     chain_atomtype, natoms)
                        dummy_atomtype_arrays = GetDummyVAtomType(atomtype_arrays)
                        lID_out.write(str(ipdb[:4])+'____chain '+str(ichain)+'____'+str(ligand_ID)+'____'+str(proteinChemName)+'____'+str(organismName)+"\n")
                    else:
                        continue

                if permutation==True:
                    rr = np.arange(natoms)
                    np.random.shuffle(rr)
                    chain_coordinate_matrix = chain_coordinate_matrix.take(rr, axis=0)
        
                print("Chain coordinate matrix shape = ",chain_coordinate_matrix.shape,
                      " and size = ",chain_coordinate_matrix.size)
                #print(chain_coordinate_matrix) # Print the Nx3 coordinate matrix of the single binding site

                dist_matrix = scipy.spatial.distance.pdist(chain_coordinate_matrix, # get the superior triangle
                                                           metric='euclidean')      # of the euclidean distance matrix
                print("Superior triangle of the distance matrix size = ", dist_matrix.size)
                dist_size_out.write(str(dist_matrix.size)+"\n")
                #print("Superior triangle of the distance matrix")
                #print(dist_matrix)
                dist_matrix_wc = np.concatenate((dist_matrix, dummy_atomtype_arrays)) 
                np.savetxt(output,dist_matrix_wc[np.newaxis],fmt='%.3f') # Save triang. dist. matrix
