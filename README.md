# AutomatedSeriesClassification

This is code for automated chemical series classification


## Original article

Automated Identification of Chemical Series: Classifying like a Medicinal Chemist
https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c00204

## Installation

The code can be downloaded and installed with:

```bash
git clone https://github.com/iwatobipen/AutomatedSeriesClassification
cd AutomatedSeriesClassification
pip install -e .
```

The ``-e`` flag means it gets installed in editable mode.

## Example usage

### Data Preparation

- run dataprep.py will download chembl27.sdf.gz and make substructurefingerprint libray.

```
$ python dataprep.py # it'll take ~30 or more minutes on my PC
```

- Then launch jupyter notebook, the notebook use same dataset as original articles. But you'll get different results compared to the article. This is because I used more newer version of ChEMBL for this code. If you would like to use same dataset to original article it is easy, just changing download link of chembl


## Acknoledgements

- Greg Landrum


## etc

Any comments, requests and suggestions will be greatly appreciated.



## License
[MIT](https://choosealicense.com/licenses/mit/)
