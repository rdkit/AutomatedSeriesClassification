#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 19:10:21 2020

@author: krugefr1
"""

from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem
import numpy as np
import matplotlib.pyplot as plt
import re

def moltosvg(mol,molSize=(450,250),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    # the MolDraw2D code is not very good at the moment at dealing with atom queries,
    # this is a workaround until that's fixed.
    # The rendering is still not going to be perfect because query bonds are not properly indicated
    opts = drawer.drawOptions()
    for atom in mc.GetAtoms():
        if atom.HasQuery() and atom.DescribeQuery().find('AtomAtomicNum')!=0:
            opts.atomLabels[atom.GetIdx()]=atom.GetSmarts()
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return svg.replace('svg:','')

def SvgsToGrid(svgs, labels, svgsPerRow=4,molSize=(250,150),fontSize=12):
    
    matcher = re.compile(r'^(<.*>\n)(<rect .*</rect>\n)(.*)</svg>',re.DOTALL) 
    hdr='' 
    ftr='</svg>' 
    rect='' 
    nRows = len(svgs)//svgsPerRow 
    if len(svgs)%svgsPerRow : nRows+=1 
    blocks = ['']*(nRows*svgsPerRow)
    labelSizeDist = fontSize*5
    fullSize=(svgsPerRow*(molSize[0]+molSize[0]/10.0),nRows*(molSize[1]+labelSizeDist))
    print(fullSize)

    count=0
    for svg,name in zip(svgs,labels):
        h,r,b = matcher.match(svg).groups()
        if not hdr: 
            hdr = h.replace("width='"+str(molSize[0])+"px'","width='%dpx'"%fullSize[0])
            hdr = hdr.replace("height='"+str(molSize[1])+"px'","height='%dpx'"%fullSize[1])
            hdr = hdr.replace("viewBox='0 0 %d %d'"%(molSize[0],molSize[1]),
                              "viewBox='0 0 %d %d'"%(fullSize[0],fullSize[1]))

        if not rect: 
            rect = r
        legend = '<text font-family="sans-serif" font-size="'+str(fontSize)+'px" text-anchor="middle" fill="black">\n'
        legend += '<tspan x="'+str(molSize[0]/2.)+'" y="'+str(molSize[1]+fontSize*2)+'">'+name.split('|')[0]+'</tspan>\n'
        if len(name.split('|')) > 1:
            legend += '<tspan x="'+str(molSize[0]/2.)+'" y="'+str(molSize[1]+fontSize*3.5)+'">'+name.split('|')[1]+'</tspan>\n'
        legend += '</text>\n'
        blocks[count] = b + legend
        count+=1

    for i,elem in enumerate(blocks): 
        row = i//svgsPerRow 
        col = i%svgsPerRow 
        elem = rect+elem 
        blocks[i] = '<g transform="translate(%d,%d)" >%s</g>'%(col*(molSize[0]+molSize[0]/10.0),row*(molSize[1]+labelSizeDist),elem) 
    res = hdr + '\n'.join(blocks)+ftr 
    return res 

def barplot_vertical(fig, ax, dict_input,ylabel,legend,xticks,legend_x_pos):
    N=len(dict_input)
    N2=len(list(dict_input.values())[0])
    ind=np.arange(N)
    width=0.35
    if N2<=20:
        clist=list(np.arange(0,20,2))+list(np.arange(1,21,2))
    else:
        clist=list(np.arange(0,N2//2*2+2,2))+list(np.arange(1,N2//2*2+1,2))
    colors=plt.cm.tab20(clist[0:len(list(dict_input.values())[0])])
    #p=np.zeros(N)
    bottom=np.zeros(N)
    for i in range(len(list(dict_input.values())[0])):
        ax.bar(ind,[v[i] for v in dict_input.values()], width, bottom=bottom, color=colors[i])
        bottom = bottom + np.array([v[i] for v in dict_input.values()])
    ax.set_ylabel(ylabel,fontsize=22)
    ax.set_xticks(ind)
    ax.set_xticklabels(xticks, rotation='vertical')
    if legend!=[]:
        ax.legend(legend, fontsize=20, ncol=4,bbox_to_anchor=(legend_x_pos, -0.3))
    plt.setp(ax.get_xticklabels(), fontsize=18)
    plt.setp(ax.get_yticklabels(), fontsize=18)
    return fig,ax