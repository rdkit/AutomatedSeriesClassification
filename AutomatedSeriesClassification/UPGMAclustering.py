#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 16:07:23 2020

@author: krugefr1
"""
import numpy as np
import sklearn
from sklearn.cluster import AgglomerativeClustering
from AutomatedSeriesClassification import utilsStructureEval
try:
    import arthor
except ImportError:
    arthor = None

def CalcSizeAndAssignment(children,Ndata):
    # Assigns molecules to the clusters of the UPGMA tree
    NumMolList=[]
    MolDict={}
    for i in range(len(children)):
        N=0
        mols_assigned=[]
        for j in range(len(children[i])):
            if children[i][j]<Ndata:
                N+=1
                mols_assigned.append(children[i][j])
            else:
                N+=NumMolList[children[i][j]-Ndata]
                mols_assigned+=MolDict[children[i][j]]
        NumMolList.append(N)
        MolDict[i+Ndata]=mols_assigned
    return NumMolList, MolDict

def CalcScore(children,distdata,NumMolList):
    # Calculates intra-cluster distance scores
    N_init=len(distdata)
    singletons=[1]*N_init
    NumMolList=singletons+NumMolList
    dist_up=np.zeros((N_init*2-1,N_init*2-1))
    dist_up[0:N_init,0:N_init]=distdata
    ScoreDict={}
    for i in range(len(children)):
        c1=children[i][0]
        c2=children[i][1]
        N1=NumMolList[c1]
        N2=NumMolList[c2]
        ScoreDict[i+N_init]=dist_up[c1,c2]
        update_x=np.zeros((1,N_init+i))
        for k in range(N_init+i):
            if k in [c1,c2]:
                dist_k=0
            else: 
                dist_k=(dist_up[c1,k]*N1+dist_up[c2,k]*N2)/(N1+N2)
            update_x[0,k]=dist_k
        dist_up[N_init+i,0:N_init+i]=update_x[0,:]
        dist_up[0:N_init+i,N_init+i]=update_x[0,:]
        
    return ScoreDict

def DetermineRelevantMCS(Ndata,children,MolDict,ScoreDict,chembldb,moldata,flimit,MinClusterSize, calcScores):
    # filter out irrelevant clusters and calculate MCS on selected clusters
    currlayer=[Ndata*2-2]
    MCSdict={}
    if hasattr(chembldb,'search'):
        Nchembl=len(chembldb.search('*'))
    else:
        Nchembl = len(chembldb)
    while len(currlayer)>0:
        childlayer=[]
        for c in currlayer:
            if c>=Ndata:
                if len(MolDict[c])>=MinClusterSize:
                    fChembl,Smarts=utilsStructureEval.MCSFromMollist(moldata.Molecule.iloc[MolDict[c]].tolist(),chembldb,Nchembl)
                    if fChembl>=flimit:
                        childlayer+=children[c-Ndata].tolist()
                    else:
                        if calcScores:
                            MCSdict[c]=(fChembl,len(MolDict[c]),Smarts,ScoreDict[c])
                        else:
                            MCSdict[c]=(fChembl,len(MolDict[c]),Smarts)
        currlayer=childlayer
    return MCSdict

def ApplyUPGMA(distdata_proj,moldata_proj,chembldb, flimit, MinClusterSize, calcScores):
    # Apply UPGMA clustering
    cluster=AgglomerativeClustering(n_clusters=2, compute_full_tree=True, affinity='precomputed',linkage='average')
    cluster.fit(distdata_proj)
    
    # Assign Clusters
    NumMolList, MolDict= CalcSizeAndAssignment(cluster.children_,len(distdata_proj))
    # Calculate intra-cluster distance scores
    if calcScores:
        ScoreDict=CalcScore(cluster.children_,distdata_proj,NumMolList)
    else: 
        ScoreDict={}
        
    # filter out irrelevant clusters and calculate MCS on selected clusters
    MCSdict=DetermineRelevantMCS(len(distdata_proj),cluster.children_,MolDict,ScoreDict,chembldb,moldata_proj,flimit,MinClusterSize, calcScores)
        
    return MCSdict


