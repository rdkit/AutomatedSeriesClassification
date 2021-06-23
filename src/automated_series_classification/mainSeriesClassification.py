#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 16:40:49 2020

@author: krugefr1
"""

import numpy as np
import os
try:
    import arthor
except ImportError:
    arthor = None
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary
import pickle
import random
import pandas as pd
import copy

from automated_series_classification import UPGMAclustering, Butinaclustering, utilsDataPrep


class Classification:
    def __init__(self,
                 proj,
                 datapath,
                 dbpath,
                 filename,
                 chembldb,
                 flimit=1e-3,
                 MinClusterSize=20,
                 clustering='UPGMA',
                 calcDists=True,
                 calcScores=False,
                 smilesCol='Smiles',
                 idCol='ID',
                 onlyCompleteRings=False,
                 useArthor=True):
        global arthor
        if not useArthor:
            arthor = None
        self.useArthor = useArthor
        self.proj = proj
        self.datapath = datapath
        self.dbpath = dbpath
        self.chembldb = chembldb
        self.flimit = flimit
        self.MinClusterSize = MinClusterSize
        self.clustering = clustering
        self.calcScores = calcScores
        self.calcDists = calcDists
        self.smilesCol = smilesCol
        self.idCol = idCol
        self.onlyCompleteRings = onlyCompleteRings
        # load data
        self.moldata_proj, self.distdata_proj = utilsDataPrep.PrepareData(
            self.proj,
            self.datapath,
            filename,
            distMeasure='Tanimoto',
            FP='Morgan2',
            calcDists=self.calcDists,
            smilesCol=smilesCol)
        if arthor is not None:
            if not os.path.isdir(dbpath):
                os.mkdir(dbpath)
            # set up project database for arthor substructure matching
            df = self.moldata_proj[[smilesCol, idCol]]
            df.to_csv('./arthor/{0}.smi'.format(self.proj),
                      header=None,
                      index=None,
                      sep=' ')
            os.system('smi2atdb -j 0 -l {0}{1}.smi {0}{1}.atdb'.format(
                self.dbpath, self.proj))
            os.system('atdb2fp -j 0 {0}{1}.atdb'.format(
                self.dbpath, self.proj))
            self.proj_db = arthor.SubDb('{0}{1}.atdb'.format(
                self.dbpath, self.proj))
        else:
            if type(dbpath) == rdSubstructLibrary.SubstructLibrary:
                self.proj_db = dbpath
                self.db_size = len(self.proj_db)
            else:
                if not os.path.exists(dbpath):
                    print("creating database")
                    mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
                    fps = rdSubstructLibrary.PatternHolder()
                    for smi in self.moldata_proj[smilesCol]:
                        m = Chem.MolFromSmiles(smi)
                        mols.AddSmiles(Chem.MolToSmiles(m))
                        fps.AddFingerprint(Chem.PatternFingerprint(m))
                    self.proj_db = rdSubstructLibrary.SubstructLibrary(
                        mols, fps)
                    self.db_size = len(mols)
                    pickle.dump(self.proj_db, open(dbpath, 'wb+'))
                else:
                    self.proj_db = pickle.load(open(dbpath, 'rb'))
                    self.db_size = len(self.proj_db)

    def AssignSeriesToMCS(self, MCSdict):
        # assign series to MCS of selected clusters
        smartslist = [v[2] for v in MCSdict.values()]
        MolAssign_prel = {}
        MolAssignment = {}
        for s in range(len(smartslist)):
            if arthor is not None:
                res = self.proj_db.search(smartslist[s])
                mols = [int(i) for i in res.to_array()]
            else:
                mols = self.proj_db.GetMatches(Chem.MolFromSmarts(
                    smartslist[s]),
                                               maxResults=self.db_size)
            MolAssign_prel[list(MCSdict.keys())[s]] = list(mols)

        # remove all series that are entirely in another series
        for key1 in MolAssign_prel.keys():
            add = 1
            for key2 in MolAssign_prel.keys():
                if key2 != key1:
                    if set(MolAssign_prel[key1]).issubset(
                            set(MolAssign_prel[key2])):
                        if set(MolAssign_prel[key2]).issubset(
                                set(MolAssign_prel[key1])) and (
                                    MCSdict[key1][0] >= MCSdict[key2][0]):
                            add = 1
                        else:
                            add = 0
                            break
            if add == 1 and MolAssign_prel[key1] not in MolAssignment.values():
                MolAssignment[key1] = MolAssign_prel[key1]

        MolAssignment = {
            k: MolAssignment[k]
            for k in MolAssignment.keys()
            if len(MolAssignment[k]) > self.MinClusterSize
        }
        if self.calcScores:
            MCSdict = {
                k: (MCSdict[k][0], len(MolAssignment[k]), MCSdict[k][2],
                    MCSdict[k][3], MolAssignment[k])
                for k in MolAssignment.keys()
            }
        else:
            MCSdict = {
                k: (MCSdict[k][0], len(MolAssignment[k]), MCSdict[k][2],
                    MolAssignment[k])
                for k in MolAssignment.keys()
            }
        return MolAssignment, MCSdict

    def ApplyClustering(self):
        # apply custering and calculate MCS
        if self.clustering == 'UPGMA':
            MCSdict = UPGMAclustering.ApplyUPGMA(
                self.distdata_proj,
                self.moldata_proj,
                self.chembldb,
                self.flimit,
                self.MinClusterSize,
                self.calcScores,
                onlyCompleteRings=self.onlyCompleteRings,
                useArthor=self.useArthor)
        elif self.clustering == 'Butina':
            distdata = copy.deepcopy(self.distdata_proj)
            MCSdict = Butinaclustering.ApplyButina(distdata,
                                                   self.moldata_proj,
                                                   self.chembldb,
                                                   self.flimit,
                                                   self.MinClusterSize,
                                                   self.calcScores,
                                                   useArthor=self.useArthor)
        else:
            print('Clustering algorithm not implemented.')
            return

        # assign series through substructure matching and filtering
        self.MolAssignment, self.MCSdict = self.AssignSeriesToMCS(MCSdict)

        # prepare and save output
        self.moldata_proj['ClusterID'] = [
            list() for x in range(self.moldata_proj.shape[0])
        ]

        for k, vs in self.MolAssignment.items():
            for v in vs:
                self.moldata_proj['ClusterID'].iloc[v].append(k)
        if self.clustering == 'UPGMA':
            self.moldata_proj.to_csv('{0}moldata_UPGMA.csv'.format(
                self.datapath))
            with open('{0}ClusterData_UPGMA.pkl'.format(self.datapath),
                      'wb') as fileout:
                pickle.dump(self.MCSdict, fileout)
        elif self.clustering == 'Butina':
            self.moldata_proj.to_csv('{0}moldata_Butina.csv'.format(
                self.datapath))
            with open('{0}ClusterData_Butina.pkl'.format(self.datapath),
                      'wb') as fileout:
                pickle.dump(self.MCSdict, fileout)
        else:
            print('Clustering algorithm not implemented.')
            return

    def CalculatePerformance(self, seriescolumn='series assignment'):

        # benchmark the automated classification against a different (probably human-defined) classification
        # human-defined compound assignment is specified in the column "seriescolumn" of the dataframe "moldata"
        # automated classification assignment specified in dict "MolAssignment"

        # calculates F1 score of automatically-identified series w.r.t. to all human-defined series, then links
        # each automatically-identified series to the human-defined series with highest F1 score

        scaflist = list(set(self.moldata_proj['scaffold'].tolist()))
        scaflist.sort()

        intersect_matrix = np.zeros((len(scaflist), len(self.MolAssignment)))
        NMatchScaf = []
        NMatchCluster = np.array([len(v) for v in self.MolAssignment.values()])
        for scaf_ind in range(len(scaflist)):
            mollist = self.moldata_proj[self.idCol].loc[self.moldata_proj[
                seriescolumn].map(lambda x: scaflist[scaf_ind] in x)].tolist()
            intersect_scaf = np.array([
                len(list(set(mollist) & set(clusterlist)))
                for clusterlist in self.MolAssignment.values()
            ])
            intersect_matrix[scaf_ind, :] = intersect_scaf
            NMatchScaf.append(len(mollist))

        NMatchScaf = np.array(NMatchScaf)
        RecallMatrix = intersect_matrix / NMatchScaf[:, None]
        PrecMatrix = intersect_matrix / NMatchCluster[None, :]
        Fscore = (2 * RecallMatrix * PrecMatrix) / (RecallMatrix + PrecMatrix +
                                                    1e-9)
        maxscore = np.argmax(Fscore, axis=0)

        PrecVector = np.zeros(len(self.MolAssignment))
        RecallVector = np.zeros(len(self.MolAssignment))
        FscoreVector = np.zeros(len(self.MolAssignment))
        LinkVector = []

        for col in range(len(self.MolAssignment)):
            PrecVector[col] = PrecMatrix[maxscore[col], col]
            RecallVector[col] = RecallMatrix[maxscore[col], col]
            FscoreVector[col] = Fscore[maxscore[col], col]
            LinkVector.append((list(self.MolAssignment.keys())[col],
                               scaflist[maxscore[col]]))

        LinkVector = np.array(LinkVector)
        self.PerformanceClusters = {
            'recall': RecallVector,
            'precision': PrecVector,
            'Fscore': FscoreVector,
            'linked series': LinkVector
        }

        if self.clustering == 'UPGMA':
            with open('{0}PerformanceData_UPGMA.pkl'.format(self.datapath),
                      'wb') as fileout:
                pickle.dump(self.PerformanceClusters, fileout)
        elif self.clustering == 'Butina':
            with open('{0}PerformanceData_Butina.pkl'.format(self.datapath),
                      'wb') as fileout:
                pickle.dump(self.PerformanceClusters, fileout)
        else:
            print('Clustering algorithm not implemented.')
            return

    def ClassificationCrossValidation(self, fraction_sample, N_sample):
        samplerange = np.arange(len(self.moldata_proj))
        invfrac = 1 / fraction_sample
        self.SampledSeries = {}
        for i in range(N_sample):

            # random sampling
            random.seed((i + 1) * 10)
            molinds = random.sample(population=samplerange.tolist(),
                                    k=int(
                                        len(samplerange.tolist()) // invfrac))
            moldata_sample = self.moldata_proj.iloc[molinds]
            distdata_sample = self.distdata_proj[molinds, :]
            distdata_sample = distdata_sample[:, molinds]

            # apply custering and calculate MCS
            if self.clustering == 'UPGMA':
                MCSdict_sampled = UPGMAclustering.ApplyUPGMA(
                    distdata_sample,
                    moldata_sample,
                    self.chembldb,
                    self.flimit,
                    self.MinClusterSize,
                    self.calcScores,
                    useArthor=self.useArthor)
            elif self.clustering == 'Butina':
                MCSdict_sampled = Butinaclustering.ApplyButina(
                    distdata_sample,
                    moldata_sample,
                    self.chembldb,
                    self.flimit,
                    self.MinClusterSize,
                    self.calcScores,
                    useArthor=self.useArthor)
            else:
                print('Clustering algorithm not implemented.')
                return

            # assign series through substructure matching and filtering
            MolAssignment_sampled, MCSdict_sampled = self.AssignSeriesToMCS(
                MCSdict_sampled)
            self.SampledSeries[i] = MCSdict_sampled

        if self.clustering == 'UPGMA':
            with open(
                    '{0}SampledSeries{1}_UPGMA.pkl'.format(
                        self.datapath, int(fraction_sample * 100)),
                    'wb') as fileout:
                pickle.dump(self.SampledSeries, fileout)
        elif self.clustering == 'Butina':
            with open(
                    '{0}SampledSeries{1}_Butina.pkl'.format(
                        self.datapath, int(fraction_sample * 100)),
                    'wb') as fileout:
                pickle.dump(self.SampledSeries, fileout)
        else:
            print('Clustering algorithm not implemented.')
            return

        return

    def EvaluationCrossValidation(self):
        # Compare the classification obtained from sampling ("SampledSeries") against the original classification ("MCSdict")
        self.EvalCrossval = pd.DataFrame(
            columns=['series id', 'repetition', 'fscore'])
        for rep in self.SampledSeries.keys():
            rep_dict = self.SampledSeries[rep]
            keylist = [k for k in rep_dict.keys()]
            for k in self.MCSdict.keys():
                intersect = [
                    len(set(self.MCSdict[k][-1]) & set(v[-1]))
                    for v in rep_dict.values()
                ]
                recall = np.array([
                    intersect[i] / len(rep_dict[keylist[i]][-1])
                    for i in range(len(keylist))
                ])
                precision = np.array(intersect) / len(self.MCSdict[k][-1])
                fscore = max(2 * recall * precision /
                             (recall + precision + 1e-9))
                row = [int(k), int(rep), fscore]
                self.EvalCrossval.loc[len(self.EvalCrossval)] = row
        self.EvalCrossval['series id'] = self.EvalCrossval['series id'].apply(
            int)
