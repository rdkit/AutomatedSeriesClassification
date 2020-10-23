#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:14:10 2020

@author: krugefr1
"""

import numpy as np
from AutomatedSeriesClassification import utilsStructureEval

def ApplyButina(distdata, moldata, chembldb, flimit, MinClusterSize,calcScores):
    MCSdict={}
    # sort assaydata ascending (i.e., lowest value is most active -> preprocess data accordingly)
    # indices: list in which molecules are selected as cluster centers
    assaydata=np.array(moldata['assay'].tolist())
    indices=np.argsort(assaydata)
    Nchembl=len(chembldb.search('*'))
    
    while len(indices)>0:
        # assign all molecules to cluster center that comply with distthresh
        # distthresh is adjusted iteratively until MCS complies with specificity threshold flimit
        distthresh=0.8
        cluster=np.where(distdata[indices[0],:]<distthresh)[0]
        if len(cluster)>=MinClusterSize:
            fChembl,Smarts=utilsStructureEval.MCSFromMollist(moldata.Molecule.iloc[cluster].tolist(),chembldb,Nchembl)
            step=0.1
            cluster_upper=cluster
            if fChembl>=flimit and len(cluster_upper)>MinClusterSize:
                while fChembl>=flimit:
                    fChembl_upper=fChembl
                    cluster_upper=cluster
                    distupper=distthresh
                    distthresh-=step
                    cluster=np.where(distdata[indices[0],:]<distthresh)[0]
                    if len(cluster)>1:
                        fChembl,Smarts=utilsStructureEval.MCSFromMollist(moldata.Molecule.iloc[cluster].tolist(),chembldb,Nchembl)
                    else: break
                distlower=distthresh
                fChembl_lower=fChembl
                fChembl_lower_old=0
                while ((fChembl_lower-fChembl_lower_old)>1e-8) and (len(cluster_upper)>=MinClusterSize):
                    distthresh=(distupper+distlower)/2
                    cluster=np.where(distdata[indices[0],:]<distthresh)[0]
                    if len(cluster)>1:
                        fChembl,Smarts=utilsStructureEval.MCSFromMollist(moldata.Molecule.iloc[cluster].tolist(),chembldb,Nchembl)
                    else: break
                    if fChembl>=flimit:
                        fChembl_upper=fChembl
                        distupper=distthresh
                        cluster_upper=cluster
                    else:
                        fChembl_lower_old=fChembl_lower
                        fChembl_lower=fChembl
                        distlower=distthresh 
            
            # select cluster/MCS if compliant with flimit and MinClusterSize
            if fChembl<flimit and len(cluster)>=MinClusterSize:
                if calcScores:
                    MCSdict[indices[0]]=(fChembl,len(cluster),Smarts,distthresh)
                else:
                    MCSdict[indices[0]]=(fChembl,len(cluster),Smarts)
        # remove molecules that were assigned to cluster (also incompliant, i.e. too small clusters)
        indices=np.array([x for x in indices if x not in cluster])
        distdata[cluster,:]=1
        distdata[:,cluster]=1
    return MCSdict