# 1. Introduction

This library is used to produce results for the official BioCreative evaluation metrics.

Written in Python 3.8

Output is printed in terminal.

# 2. Requirements

+ Python3
+ pandas

You will need python3 (together with its base libraries) and the pandas package version.

To install them: 

```
git clone URL
cd drugprot-evaluation-library
pip install -r requirements.txt
```


# 3. Execution

```
cd src
python main.py -g ../gs-data/gs_relations.tsv -p ../toy-data/pred_relations.tsv -e ../gs-data/gs_entities.tsv --pmids ../gs-data/pmids.txt
```

To run the evaluation library, move to the src/ directory and execute the main.py script.


# 4. Other interesting stuff:
### Metrics
The relevant metrics are micro-average precision, recall and f1-score. The tool allows you to explore your results by relation type as well.


### Data Format
The **Predictions TSV file contains the predicted relations**. It *must* have these columns (separated by a \t):
+ Article identifier (PMID)
+ DrugProt relation
+ Interactor argument 1 (of type CHEMICAL)
+ Interactor argument 2 (of type GENE)

Example:

12488248	  &nbsp;&nbsp; INHIBITOR	  &nbsp;&nbsp; Arg1:T1	  &nbsp;&nbsp; Arg2:T52

12488248	  &nbsp;&nbsp; INHIBITOR	  &nbsp;&nbsp; Arg1:T2	  &nbsp;&nbsp; Arg2:T52

23220562	  &nbsp;&nbsp; ACTIVATOR	  &nbsp;&nbsp; Arg1:T12	  &nbsp;&nbsp; Arg2:T42

23220562	  &nbsp;&nbsp; ACTIVATOR	  &nbsp;&nbsp; Arg1:T12	  &nbsp;&nbsp; Arg2:T43

23220562	  &nbsp;&nbsp; INDIRECT-DOWNREGULATOR	  &nbsp;&nbsp; Arg1:T1	  &nbsp;&nbsp; Arg2:T14


For more in-depth information about the Data Format (Gold Standard and Predictions), have a look at the [toy-data](toy-data) directory or at the [Zenodo page](https://doi.org/10.5281/zenodo.4955410) .

### Script Arguments
+ ```-g/--gs_path```: path to Gold Standard relations TSV file
+ ```-p/--pred_path```: path to Prediction TSV file
+ ```-e/--ent_path```: path to Gold Standard entities TSV file
+ ```--pmids```: path to list of relevant PMIDs


### Examples: 

```
$ cd src
$ python main.py -g ../gs-data/gs_relations.tsv -p ../toy-data/pred_relations.tsv -e ../gs-data/gs_entities.tsv --pmids ../gs-data/pmids.txt

python main.py -g ../gs-data/gs_relations.tsv -p ../toy-data/pred_relations.tsv -e ../gs-data/gs_entities.tsv --pmids ../gs-data/pmids.txt
Loading GS files...
Loading prediction files...
Checking GS files...
Checking Predictions files...
Formatting data...
Computing DrugProt (BioCreative VII) metrics ...
(p = Precision, r=Recall, f1 = F1 score)
By relation type
p_INDIRECT-DOWNREGULATOR=1.0
r_INDIRECT-DOWNREGULATOR=0.5
f1_INDIRECT-DOWNREGULATOR=0.667
p_DIRECT-REGULATOR=1.0
r_DIRECT-REGULATOR=1.0
f1_DIRECT-REGULATOR=1.0
p_ACTIVATOR=1.0
r_ACTIVATOR=1.0
f1_ACTIVATOR=1.0
p_INHIBITOR=1.0
r_INHIBITOR=1.0
f1_INHIBITOR=1.0
p_AGONIST=0.8
r_AGONIST=0.667
f1_AGONIST=0.727
p_ANTAGONIST=1.0
r_ANTAGONIST=1.0
f1_ANTAGONIST=1.0
p_PART-OF=1.0
r_PART-OF=0.5
f1_PART-OF=0.667
The following relations are not present in the Gold Standard: INDIRECT-UPREGULATOR,PRODUCT-OF,AGONIST-INHIBITOR,SUBSTRATE_PRODUCT-OF,AGONIST-ACTIVATOR,SUBSTRATE
The following relations are not present in the Predictions: SUBSTRATE_PRODUCT-OF,PRODUCT-OF,AGONIST-ACTIVATOR,SUBSTRATE

Gobal results across all DrugProt relations (micro-average)
p_micro=0.783
r_micro=0.75
f1_micro=0.766

```

### Relevant links:
+ [DrugProt Web](https://biocreative.bioinformatics.udel.edu/tasks/biocreative-vii/track-1/)
+ [DrugProt corpus](https://doi.org/10.5281/zenodo.4955410)
