import gzip
import os
import pickle

import click
import wget
from rdkit import Chem, RDLogger, rdBase
from tqdm import tqdm

RDLogger.DisableLog("rdApp.warning")

# See https://pubs.acs.org/doi/10.1021/jm020472j
bradley_url = 'https://pubs.acs.org/doi/suppl/10.1021/jm020472j/suppl_file/jm020472j_s2.xls'


@click.command()
@click.option(
    '--directory', type=click.Path(file_okay=False, dir_okay=True),
    default=os.getcwd, help='Defaults to current directory',
)
@click.option('--chebml-version', default='27', show_default=True)
def main(directory: str, chebml_version: str):
    """Download the ChEBML data."""
    os.makedirs(directory, exist_ok=True)

    bradley_path = os.path.join(directory, 'jm020472j_s2.xls')
    if not os.path.exists(bradley_path):
        try:
            wget.download(bradley_url, out=directory)
        except:
            click.echo('There goes ACS stopping science')

    chembl_url = (
        f'ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/'
        f'chembl_{chebml_version}/chembl_{chebml_version}.sdf.gz'
    )

    sdf_path = os.path.join(directory, f'chembl_{chebml_version}.sdf.gz')
    if not os.path.exists(sdf_path):
        wget.download(chembl_url, out=directory)

    sss_path = os.path.join(directory, f'chembl{chebml_version}_sssdata.pkl')
    if not os.path.exists(sss_path):
        click.echo(f'RDKit Version: {rdBase.rdkitVersion}')
        data = []

        with gzip.GzipFile(sdf_path) as gz:
            suppl = Chem.ForwardSDMolSupplier(gz)
            for mol in tqdm(suppl, desc=f'Processing ChEBML {chebml_version}', unit_scale=True):
                if mol is None or mol.GetNumAtoms() > 50:
                    continue
                fp = Chem.PatternFingerprint(mol)
                smi = Chem.MolToSmiles(mol)
                data.append((smi, fp))

        click.echo(f'Outputting to {sss_path}')
        with open(sss_path, 'wb') as file:
            pickle.dump(data, file, protocol=pickle.HIGHEST_PROTOCOL)

    click.echo('Done ;)')


if __name__ == '__main__':
    main()
