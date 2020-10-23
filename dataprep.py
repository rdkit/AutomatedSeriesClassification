import gzip
import os
import pickle

import wget
from rdkit import Chem, RDLogger, rdBase
from tqdm import tqdm

print('RDKit Version', rdBase.rdkitVersion)
RDLogger.DisableLog("rdApp.warning")

url = 'ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_27.sdf.gz'
sdf_path = './chembl_27.sdf.gz'
sss_path = './data/chembl27_sssdata.pkl'

if not os.path.exists(sdf_path):
    filename = wget.download(url)

os.makedirs('./data', exist_ok=True)

if not os.path.exists(sss_path):
    data = []

    with gzip.GzipFile(sdf_path) as gz:
        suppl = Chem.ForwardSDMolSupplier(gz)
        for mol in tqdm(suppl, desc='Processing molecules', unit_scale=True):
            if mol is None or mol.GetNumAtoms() > 50:
                continue
            fp = Chem.PatternFingerprint(mol)
            smi = Chem.MolToSmiles(mol)
            data.append((smi, fp))

    with open(sss_path, 'wb') as file:
        pickle.dump(data, file, protocol=pickle.HIGHEST_PROTOCOL)

print('Done ;)')
