import os
import wget
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary
from rdkit import RDLogger
from rdkit import rdBase
import pickle
import time
import gzip
print(rdBase.rdkitVersion)
RDLogger.DisableLog("rdApp.warning")

url = 'ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_27.sdf.gz'
if not os.path.exists('./chembl_27.sdf.gz'):
    filename = wget.download(url)

gz = gzip.GzipFile('./chembl_27.sdf.gz')
suppl = Chem.ForwardSDMolSupplier(gz)


if not os.path.exists('./data'):
    os.mkdir('./data')
if not os.path.exists('./data/chembl27_sssdata.pkl'):
    t1=time.time()
    data = []
    for i,mol in enumerate(suppl):
        if not ((i+1)%50000):
            print(f"Processed {i+1} molecules in {(time.time()-t1):.1f} seconds")
        if mol is None or mol.GetNumAtoms()>50:
            continue
        fp = Chem.PatternFingerprint(mol)
        smi = Chem.MolToSmiles(mol)
        data.append((smi,fp))
    t2=time.time()
    pickle.dump(data,open('./data/chembl27_sssdata.pkl','wb+'))

print('Done ;)')
