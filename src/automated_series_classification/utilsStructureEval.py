#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 16:29:23 2020

@author: krugefr1
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdSubstructLibrary
import pickle
try:
    import arthor
except ImportError:
    arthor = None


def MCSFromMollist(mollist,chembldb,Nchembl):
    MCSSmarts2=rdFMCS.FindMCS(mollist,atomCompare=rdFMCS.AtomCompare.CompareAny,bondCompare=rdFMCS.BondCompare.CompareOrderExact,ringMatchesRingOnly=False,timeout=1).smartsString
    MCSSmarts=rdFMCS.FindMCS(mollist,atomCompare=rdFMCS.AtomCompare.CompareElements,bondCompare=rdFMCS.BondCompare.CompareOrder,ringMatchesRingOnly=False,timeout=1).smartsString
    if MCSSmarts2=='': fChembl2=1
    else: fChembl2=getFChembl(MCSSmarts2,chembldb,Nchembl)
    if MCSSmarts=='': fChembl=1
    else:fChembl=getFChembl(MCSSmarts,chembldb,Nchembl)
    if fChembl2<fChembl:
        fChembl=fChembl2
        MCSSmarts=MCSSmarts2
    return fChembl,MCSSmarts


def _arthor_getFChembl(qry,chembldb,Ntot,qryformat='Smarts'):
    if qryformat=='Smarts':
        results=chembldb.search(qry)
    elif qryformat=='MDL':
        with open(qry) as f:
            qryarthor=arthor.Query(f.read(),"Mdl")
        results=chembldb.search(str(qryarthor))
    fChembl=(len(results)+1)/(Ntot+2)
    return fChembl 


_ssslib=None
def _rdkit_getFChembl(qry,chembldb,Ntot,qryformat='Smarts'):
    global _ssslib
    if chembldb is None:
        if _ssslib is None:
            _ssslib = pickle.load(open('./chembl27_sslib.pkl','rb'))
        chembldb = _ssslib
    if qryformat=='Smarts':
        qry = Chem.MolFromSmarts(qry)
    elif qryformat=='MDL':
        qry = Chem.MolFromMolFile(qry)
    results=chembldb.GetMatches(qry,maxResults=Ntot)

    fChembl=(len(results)+1)/(Ntot+2)
    return fChembl 


def getFChembl(qry,chembldb,Ntot,qryformat='Smarts'):
    if arthor is not None:
        return _arthor_getFChembl(qry,chembldb,Ntot,qryformat=qryformat)
    else:
        return _rdkit_getFChembl(qry,chembldb,Ntot,qryformat=qryformat)