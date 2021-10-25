#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 15:32:08 2021

@author: tonifuc3m
"""

import pandas as pd
import warnings
from sklearn.metrics import f1_score, precision_score, recall_score
import numpy as np
# from numpy import mean

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
warnings.formatwarning = warning_on_one_line

def Flatten(ul):
    fl = []
    for i in ul:
        if type(i) is list:
            fl += Flatten(i)
        else:
            fl += [i]
    return fl

def main(y_true, y_pred, reltype2tag, gs_rel_list, pred_rel_list, pmid2pos):
    '''
    Compute precision, recall and F1-score and print them.

    Parameters
    ----------
    y_true : list
        List of GS Relations.
    y_pred : list
        List of Predictions Relations.
    reltype2tag : dict
        Mapping from relation type string to integer tag
    gs_rel_list : list
        List of relation types in GS
    pred_rel_list : list
        List of relation types in Predictions
    pmid2pos : dict
        Dictionary with PMIDs and they positions in y_true and y_pred of the
        annotations corresponding to each PMID

    Returns
    -------
    None.

    '''
    relations_not_in_gs = set(reltype2tag.keys()) - set(gs_rel_list)
    relations_not_in_pred = set(reltype2tag.keys()) - set(pred_rel_list)
    
    
    print("Number of predictions per relation type")
    for k,v in sorted(reltype2tag.items()):          
        y_pred_this = y_pred[v-1]
        print(f"{k}={sum(y_pred_this)}")
        
    print("Number of GS relations per relation type")
    for k,v in sorted(reltype2tag.items()):          
        y_true_this = y_true[v-1]
        print(f"{k}={sum(y_true_this)}")
        
        
    print("By relation type")
    f1_by_rel_type = []
    p_by_rel_type = []
    r_by_rel_type = []
    for k,v in sorted(reltype2tag.items()):
        if k in relations_not_in_gs:
            continue
        y_true_this = y_true[v-1]
        y_pred_this = y_pred[v-1]
        f1_this = f1_score(y_true_this, y_pred_this, zero_division=0)
        p_this = precision_score(y_true_this, y_pred_this, zero_division=0)
        r_this = recall_score(y_true_this, y_pred_this, zero_division=0)
        f1_by_rel_type.append(f1_this)
        p_by_rel_type.append(p_this)
        r_by_rel_type.append(r_this)
        print(f"p_micro_{k}={round(p_this, 6)}\nr_micro_{k}={round(r_this, 6)}\nf1_micro_{k}={round(f1_this, 6)}")
        
    print(f"The following relations are not present in the Gold Standard: {','.join(relations_not_in_gs)}")
    print(f"The following relations are not present in the Predictions: {','.join(relations_not_in_pred)}")

    print("\nGobal results across all DrugProt relations (micro-average)")
    F1 = f1_score(y_true, y_pred, average='micro')
    P = precision_score(y_true, y_pred, average='micro')
    R = recall_score(y_true, y_pred, average='micro')
    print(f"p_micro={round(P, 3)}\nr_micro={round(R, 3)}\nf1_micro={round(F1, 3)}")
    
    print("\nMacro-average (per PMID)")
    f1_per_pmid = []
    p_per_pmid = []
    r_per_pmid = []
    for pmid, pos in pmid2pos.items():
        if pos[0]==pos[1]:
            continue
        y_true_this = list(map(lambda x: x[pos[0]:(pos[1])], y_true))
        y_pred_this = list(map(lambda x: x[pos[0]:(pos[1])], y_pred))
        f1_per_pmid.append(f1_score(y_true_this, y_pred_this, average='micro'))
        p_per_pmid.append(precision_score(y_true_this, y_pred_this, average='micro'))
        r_per_pmid.append(recall_score(y_true_this, y_pred_this, average='micro'))
        # for k,v in reltype2tag.items():
        #     print(f"GS Positives {k}= {sum(Flatten(y_true_this[v-1]))}")
        #     print(f"Pred Positives {k}= {sum(Flatten(y_pred_this[v-1]))}")
        # print(p_per_pmid, r_per_pmid, f1_per_pmid)

    print(f"p_macro={round(np.mean(p_per_pmid), 3)}\nr_macro={round(np.mean(r_per_pmid), 3)}\nf1_macro={round(np.mean(f1_per_pmid), 3)}")
