# Class BindingDB	
# v. 0.0
# G. Tiana, 18 May 2020
#

import pandas as pd
import sys

class BindingDB:

   def __init__(self, filename, verbose=True, dropRepeated=False, flattenPdb=False):

     if verbose:
       print("Reading entries from "+filename)

     try: 
       self.df = pd.read_csv(filename, sep='\t', header=0,  dtype=str )
     except:
       sys.exit("ERROR: File "+filename+" not found.")


     if dropRepeated:
       self.df.drop_duplicates( subset=['PubChem SID','UniProt (SwissProt) Primary ID of Target Chain'], keep='first', inplace=True )


     # split multiple pdbs
     self.df['PDB ID(s) for Ligand-Target Complex'] = self.df['PDB ID(s) for Ligand-Target Complex'].apply( lambda x: str(x).split(',')  ) 

     if flattenPdb:
       self.df = self.df.explode('PDB ID(s) for Ligand-Target Complex', ignore_index=True)

     if verbose:
       print(" "+str(len(self.df.index))+" entries read.") 

