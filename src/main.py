#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 15:32:08 2021

@author: tonifuc3m
"""

import argparse
import warnings
import pandas as pd
import os

import compute_metrics
from utils import load_entities_dict, prepro_relations, \
    format_relations, get_chemical_gene_combinations

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
warnings.formatwarning = warning_on_one_line

def parse_arguments():
    '''
    DESCRIPTION: Parse command line arguments
    '''
  
    parser = argparse.ArgumentParser(description='process user given parameters')
    parser.add_argument('-g', '--gs_path', required = False, dest = 'gs_path', 
                        default = '../gs-data/gs_relations.tsv', 
                        help = 'path to GS relations file (TSV)')
    parser.add_argument('-e', '--ent_path', required = False, dest = 'ent_path', 
                        default = '../gs-data/gs_entities.tsv', 
                        help = 'path to GS entities file (TSV)')
    parser.add_argument('-p', '--pred_path', required = False, dest = 'pred_path', 
                        default = '../toy-data/pred_relations.tsv', 
                        help = 'path to predictions file (TSV)')
    parser.add_argument('--pmids', required = False, dest = 'pmids', 
                        default = '../gs-data/pmids.txt',
                        help = 'path to list of valid pubmed IDs. One PMID per line')
    
    return parser.parse_args()

def main(args):
    '''
    Load GS and Predictions; format them; compute precision, recall and 
    F1-score and print them.

    Parameters
    ----------
    gs_path : str
        Path to GS Relations TSV file.
    pred_path : str
        Path to Predictions Relations TSV file.
    ent_path : str
        Path to GS Entities TSV file
    pmids : str
        Path to file with valid pubmed IDs

    Returns
    -------
    None.

    '''
    
    rel_types = ['INDIRECT-DOWNREGULATOR','INDIRECT-UPREGULATOR','DIRECT-REGULATOR',
             'ACTIVATOR','INHIBITOR','AGONIST','AGONIST-ACTIVATOR',
             'AGONIST-INHIBITOR','ANTAGONIST','PRODUCT-OF','SUBSTRATE',
             'SUBSTRATE_PRODUCT-OF','PART-OF']
    reltype2tag = {w: i+1 for i, w in enumerate(rel_types)}
    NREL = len(reltype2tag.keys())
        
    # Load GS
    print("Loading GS files...")
    _dict_, genes, chemicals = load_entities_dict(args.ent_path)
    combinations, NCOMB = get_chemical_gene_combinations(_dict_)
    pmids = set(map(lambda x: str(x.strip()), open(args.pmids)))
    gs = pd.read_csv(args.gs_path, sep='\t', header=None, dtype=str, skip_blank_lines=True,
                     names = ['pmid', 'rel_type', 'arg1', 'arg2'], encoding = 'utf-8')
    
    # Load predictions
    print("Loading prediction files...")
    pred = pd.read_csv(args.pred_path, sep='\t', header=None, dtype=str, skip_blank_lines=True,
                       names = ['pmid', 'rel_type', 'arg1', 'arg2'], encoding = 'utf-8')
    
    # Format data
    print("Checking GS files...")
    gs_valid,gs_rel_list = prepro_relations(gs, chemicals, rel_types, is_gs=True)
    
    print("Checking Predictions files...")
    pred_valid,pred_rel_list = prepro_relations(pred, chemicals, rel_types, is_gs=False, gs_files=pmids)
    
    print("Formatting data...")
    y_true, y_pred = format_relations(gs_valid, pred_valid, combinations, 
                                      NCOMB, NREL, reltype2tag)
    
    # Compute metrics
    print("Computing DrugProt (BioCreative VII) metrics ...\n(p = Precision, r=Recall, f1 = F1 score)")
    compute_metrics.main(y_true, y_pred, reltype2tag, gs_rel_list, pred_rel_list)
    
if __name__ == '__main__':
    
    args = parse_arguments()
    
    if os.path.exists(args.gs_path)==False:
        raise Exception(f'Gold Standard path {args.gs_path} does not exist')
    if os.path.exists(args.pred_path)==False:
        raise Exception(f'Predictions path {args.pred_path} does not exist')
    if os.path.exists(args.ent_path)==False:
        raise Exception(f'Gold Standard entities path {args.ent_path} does not exist')
    if os.path.exists(args.pmids)==False:
        raise Exception(f'PMIDs file list path {args.pmids} does not exist')
        
    main(args)
