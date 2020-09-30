#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 14:51:39 2020

@author: krugefr1
"""

from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools, DataStructs
import pickle
import pandas as pd
import numpy as np

def calcDistMatrix(df, distMeasure):
    # calculates the distance matrix between all paris of molecules, standard: Tanimoto and Morgan2 FPs
    dists=np.zeros([len(df),len(df)])
    if distMeasure=='Tanimoto':
        for i in range(1,len(df)):
            ds = DataStructs.BulkTanimotoSimilarity(df.FP.iloc[i],list(df.FP.iloc[:i]),returnDistance=1)
            for j in range(i):
                dists[i,j] = ds[j]
                dists[j,i] = ds[j]
    else:
        print(distMeasure, 'distance metric not implemented.')
        return
    return dists 

def readProjectData(filename, FP):
    # reads in the project data and calculates fingerprints
    df_proj=pd.read_csv(filename,names=['ID','Structure','mol name','scaffold','series assignment','assay'], skiprows=[0])
    #df_proj = df_proj.head(100)
    PandasTools.AddMoleculeColumnToFrame(df_proj,smilesCol='Structure',molCol='Molecule')
    df_proj=df_proj.loc[df_proj['Molecule'].map(lambda x: x is not None)]
    if FP=='Morgan2':
        df_proj['FP']=df_proj.Molecule.map(lambda x : AllChem.GetMorganFingerprint(x,2))
    else: 
        print(FP, ' fingerprint not implemented.')
        return
    return df_proj


def PrepareData(proj,datapath,distMeasure='Tanimoto',FP='Morgan2', calcDists=False):
    # reads in project data and distance matrix (or calculate distance matrix)
    filename='{0}moldata_preprocessed.csv'.format(datapath)
    moldata=readProjectData(filename, FP)
    print(f'read {len(moldata)} molecules')
    if calcDists:
        dists=calcDistMatrix(moldata, distMeasure)
        with open('{0}distmatrix.txt'.format(datapath), 'wb') as fileout:
            pickle.dump(dists,fileout)
    else:
        with open('{0}distmatrix.txt'.format(datapath), 'rb') as filein:
            dists=pickle.load(filein)
    
    return moldata, dists
