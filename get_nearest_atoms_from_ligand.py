#!/usr/bin/env python3.7

from PDB import Pdb
from bindingDB import BindingDB
import requests
import sys
from os import path,mkdir

# Functions
def DownloadPdb(pdbCode, dirPdb='.', verbose=True):

  if verbose: print('Downloading '+pdbCode)

  url = 'https://files.rcsb.org/download/'+pdbCode+'.pdb'
  try:
    r = requests.get(url, allow_redirects=True)
    open(dirPdb+'/'+pdbCode+'.pdb', 'wb').write(r.content)
  except:
    sys.exit('ERROR: cannot download pdb file')

  if verbose: print(' saved to '+dirPdb+'/')
################################################################################

dirPdb = 'pdbs'
dirNewPdb = 'nearest_chains_pdbs'
verbose = True
distanceContact = 5     #in Angstrom

pdbs_stat = open("pdbs_statistics.dat","w")

# Define database
db = BindingDB("PDB_BindingDB.tsv", dropRepeated=True, flattenPdb=True)
# keep rows that has no NaN value in these columns
db.df = db.df[db.df['PDB ID(s) for Ligand-Target Complex'].notna()] 
db.df = db.df[db.df['Ligand HET ID in PDB'].notna()] 
# Removing rows with 'nan' string as pdbcode
db.df = db.df[~db.df['PDB ID(s) for Ligand-Target Complex'].str.contains('nan')]
db.df = db.df[~db.df['Ligand HET ID in PDB'].str.contains('nan')]

pdbList = db.df['PDB ID(s) for Ligand-Target Complex']
ligandHETID = db.df['Ligand HET ID in PDB']
# protein name
proteinChemName =db.df['UniProt (SwissProt) Recommended Name of Target Chain']
proteinChemName.fillna("unknown", inplace = True) 
# organism name
organismName = db.df['UniProt (SwissProt) Entry Name of Target Chain']
organismName = organismName.str.split('_').str[-1]
organismName.fillna("unknown", inplace = True) 

# Check dirs
if not path.exists(dirPdb):
  try:
    mkdir(dirPdb)
  except:
    sys.exit('ERROR: Cannot create directory for pdb')

if not path.exists(dirNewPdb):
    try:
        mkdir(dirNewPdb)
    except:
        sys.exit('ERROR: Cannot create directory for pdb of nearest chains')

# Loop on list of pdb
for line1, line2, line3, line4 in zip(pdbList, ligandHETID, proteinChemName, organismName):

    pdbCode = line1.strip()
    ligandHETID = line2.strip()
    proteinChemName = line3.strip()
    organismName = line4.strip()
    if verbose: print('Analyzing: '+str(pdbCode)+', ligand HET ID: '+str(ligandHETID)+', protein chem-name: '+str(proteinChemName)+', organism: '+str(organismName))

    # download if not there
    fileName = dirPdb+'/'+pdbCode+'.pdb'
    if not path.isfile(fileName):
        DownloadPdb(pdbCode,dirPdb)

    pdb = Pdb()
    pdb.verbose = True
    pdb.ReadFile(fileName,ligandHETID)

    # loop on disjoint chains
    chainString = pdb.chains
    protString = ''
    chemString = ''

    for c in chainString:
        chainType = pdb.IdentifyChainType(c)

        if chainType == 'chemical':
            chemString = chemString + c

        if chainType == 'protein':
            protString = protString + c

    if len(chemString)>1:
        print("WARNING: more than one chemical chain ("+chemString+")")

    if pdb.ligandHETID !='':
        pdb.RemoveNotHETIDAtomsFromChemicalChains(chemString)
    else:
        print("WARNING: pdb has no HET ID for the ligand. Discarding pdb.")
        continue

    #for c in chainString:
    #    chainType = pdb.IdentifyChainType(c)
    
    # select all the atoms with dist<= distance
    selPdb = pdb.SelectChainByDistance(chainFromList=protString, chainRefList=chemString, distance=distanceContact)
    selPdb.Print(dirNewPdb+'/'+pdbCode+'_nearest_chains_to_ligand.pdb')

    pdbs_stat.write(str(pdbCode)+'____'+str(ligandHETID)+'____'+str(proteinChemName)+'____'+str(organismName)
                    +'____TOTAL_NUMBER_ATOMS____'+str(pdb.nAtoms)+'____SELECTED_ATOMS____'
                    +str(selPdb.nAtoms)+"\n")

    if verbose: print()
    del pdb
    del selPdb

pdbs_stat.close()
