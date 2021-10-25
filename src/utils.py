#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 14:32:25 2021

@author: antonio
"""
import warnings
import itertools


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
warnings.formatwarning = warning_on_one_line

def save_mark(chemicals, genes, ent_type, mark):
    if ent_type=='CHEMICAL':
        chemicals.append(mark)
    elif ent_type=='GENE':
        genes.append(mark)
    else:
        warnings.warn('Wrong Entity type')
    return chemicals, genes

def update_dict(chemicals, genes, pmid, _dict_):
    _dict_this = {}
    _dict_this['chemicals'] = chemicals
    _dict_this['genes'] = genes
    _dict_[pmid] = _dict_this
    return _dict_

def load_entities_dict(path):
    """
    Load entities TSV
    
    Returns
    -------
    genes : set
        Set of GENE
    chemicals : set
        Set of CHEMICAL
    _dict_ : dict
        Dictionary with annotated entities. Keys: PMID, Values: another
        dictionary with keys 'chemicals' and 'genes' and value the 
        annotation mark
    """
    _dict_ = {}
    general_dict = {}
    pmid_prev = ''
    chemicals = []
    genes = []
    for line in open(path).readlines():
        split = line.split('\t')
        if (len(line.split('\t')) != 6) & (line!='\n'):
            raise Exception(f"Line {line} in file {path} wrongly formatted")
        pmid = split[0]
        mark = split[1]
        ent_type = split[2]
        
        if pmid != pmid_prev:
            # Save previous PMID
            _dict_ = update_dict(chemicals, genes, pmid_prev, _dict_)
            
            # Initialize new PMID
            chemicals = []
            genes = []
        
        # Save entity
        chemicals, genes = save_mark(chemicals, genes, ent_type, mark)
        
        # Update previous PMID
        pmid_prev = pmid
        
        # Save General dict
        general_dict[pmid + '-' + mark] = ent_type
        
    # Save last PMID
    _dict_ = update_dict(chemicals, genes, pmid_prev, _dict_)
    del _dict_['']
    
    # Get list of genes and chemicals
    genes = set([k for k,v in general_dict.items() if v=='GENE'])
    chemicals = set([k for k,v in general_dict.items() if v=='CHEMICAL'])
    
    return _dict_, genes, chemicals

def get_chemical_gene_combinations(_dict_):
    """
    Parameters
    ----------
    _dict_ : dictionary
        Dictionary with annotated entities. Keys: PMID, Values: another
        dictionary with keys 'chemicals' and 'genes' and value the 
        annotation mark

    Returns
    -------
    combinations : dictionary
        PMIDs as keys, all possible CHEMICAL-GENE combinations are values.
    NCOMB : int
        DESCRIPTION.

    """

    combinations = {}
    pmid2pos = {}
    NCOMB = 0
    pos1 = 0
    for pmid, entities in _dict_.items():
        chem = entities['chemicals']
        genes = entities['genes']
        combinations[pmid] = list(itertools.product(chem, genes))
        NCOMB = NCOMB + len(combinations[pmid] )
        pos0 = pos1
        pos1 = pos1 + len(combinations[pmid])
        pmid2pos[pmid] = (pos0, pos1)
    return combinations, NCOMB, pmid2pos

def prepro_relations(df, chemicals, rel_types, is_gs=False, gs_files=set()):
    """
    Preprocess annotations dataframe

    Parameters
    ----------
    df : pandas DataFrame
        Dataframe with annotations (GS or predicted).
    chemicals : list
        List of GS CHEMICAL entities.
    rel_types : list
        List of valid relation types.
    is_gs : bool, optional
        Whether we are formatting the GS annotations. The default is False.
    gs_files : set, optional
        Set of PMIDs. The default is set().

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    df : pandas DataFrame
        Clean annotations DataFrame.

    """
    
    if df.shape[0] == 0:
        raise Exception('There are not parsed annotations')
    if df.shape[1] != 4:
        raise Exception('Wrong column number in the annotations file')
            

    # Remove predictions for PMIDs not valid
    if is_gs==False:
        df = df.loc[df['pmid'].isin(gs_files),:].copy()
        
    # Remove predictions for RELATION FILES not valid
    if len(set(df.rel_type.tolist()).intersection(set(rel_types))) > len(set(rel_types)):
        warnings.warn("Non-valid relation types. Skipping them")
    df = df.loc[df['rel_type'].isin(rel_types),:].copy()
    
    # Drop duplicates
    df = df.drop_duplicates(subset=df.columns).copy()
    
    # Check every relation has one CHEMICAL and one GENE
    df['pmid-arg1'] = df['pmid'] + '-' + df['arg1'].apply(lambda x: x.split(':')[-1])
    df['pmid-arg2'] = df['pmid'] + '-' + df['arg2'].apply(lambda x: x.split(':')[-1])
    df['is_arg1_chemical'] = df['pmid-arg1'].apply(lambda x: x in chemicals)
    df['is_arg2_chemical'] = df['pmid-arg2'].apply(lambda x: x in chemicals)
    skip_chem = []
    skip_gene = []
    if any(df.apply(lambda x: x['is_arg1_chemical'] + x['is_arg2_chemical'], axis=1) > 1):
        skip_chem = df.loc[df.apply(lambda x: x['is_arg1_chemical'] + x['is_arg2_chemical'],
                               axis=1) > 1].index.tolist()
        warnings.warn(f"The following lines have more than one CHEMICAL entity: {df.loc[skip_chem]}. Skipping them")
    if any(df.apply(lambda x: x['is_arg1_chemical'] + x['is_arg2_chemical'], axis=1) == 0):
        skip_gene = df.loc[df.apply(lambda x: x['is_arg1_chemical'] + x['is_arg2_chemical'],
                               axis=1) == 0].index.tolist()
        warnings.warn(f"The following lines have less than one CHEMICAL entity: {df.loc[skip_gene]}. Skipping them")
    skip = skip_chem + skip_gene
    if len(skip)>1:
        df.drop(skip, inplace=True)


    # Get final CHEMICAL and gene
    # TODO: double-check this step has no errors
    df['chemical'] = df.apply(lambda x: x['arg1'].split(':')[-1] \
                              if x['is_arg1_chemical']==True else x['arg2'].split(':')[-1], axis=1)
    df['gene'] = df.apply(lambda x: x['arg2'].split(':')[-1] \
                              if x['is_arg1_chemical']==True else x['arg1'].split(':')[-1], axis=1)
        
    # Keep only relevant columns
    df = df[['pmid', 'rel_type', 'chemical', 'gene']].drop_duplicates(subset=['pmid', 'rel_type', 'chemical', 'gene']).copy()
    
    return df, set(df.rel_type.tolist())


def format_relations(gs_valid, pred_valid, combinations, NCOMB, NREL, reltype2tag):
    """
    Format relation information for sklearn

    Parameters
    ----------
    gs_valid : pandas DataFrame
        GS relations. Columns: ['pmid', 'rel_type', 'chemical', 'gene']
    pred_valid : pandas DataFrame
        Predicted relations. ['pmid', 'rel_type', 'chemical', 'gene']
    combinations : dictionary
        PMIDs as keys, all possible CHEMICAL-GENE combinations are values.
    NCOMB : int
        Number of entity combinations (CHEMICAL-GENE).
    NREL : int
        Number of valid relations.
    reltype2tag : dict
        Mapping from relation type string to integer tag.

    Returns
    -------
    y_true : nested lists
        List of GS relations.
    y_pred : nested lists
        List of Predictions relations.

    """
    
    y_true = [[0]*NCOMB for _ in range(NREL)]
    y_pred = [[0]*NCOMB for _ in range(NREL)]
    pos0 = 0
    for pmid,combs in sorted(combinations.items()):
        if combs==[]:
            continue
        
        # Subset GS and predictions
        gs_this = gs_valid.loc[gs_valid['pmid']==pmid,:]
        pred_this = pred_valid.loc[pred_valid['pmid']==pmid,:]
        
        # Iterate over all combinations
        for c, idx in zip(combs, range(len(combs))):
            chem = c[0]
            gene = c[1]
            
            gs_rel = gs_this.loc[(gs_this['chemical']==chem)&(gs_this['gene']==gene),
                                   'rel_type'].values
            if len(gs_rel)>0:
                for item in gs_rel:
                    tag = reltype2tag[item]
                    y_true[tag-1][pos0+idx]=1
                
            pred_rel = pred_this.loc[(pred_this['chemical']==chem)&(pred_this['gene']==gene),
                                   'rel_type'].values
            if len(pred_rel)>0:
                for item in pred_rel:
                    tag = reltype2tag[item]
                    y_pred[tag-1][pos0+idx]=1
        pos0 = pos0+idx+1
    return y_true, y_pred


def filter_pred(pred_path, pmids):
    '''
    Create temporary predictions file with only the PMIDs that are relevant
    for evaluation

    Parameters
    ----------
    pred_path : str
        Predictions file.
    pmids : set
        Set of PMIDs relevant for evaluation.

    Returns
    -------
    outpath : str
        Temporary predictions file.

    '''
    outpath = pred_path + '.tmp.tsv'
    fout = open(outpath, 'w', encoding='utf-8')
    with open(pred_path) as fin:
        for line in fin:
            if line.split('\t')[0] in pmids:
                fout.write(line)
                
        
    return outpath