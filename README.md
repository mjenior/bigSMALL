BiGSMAll or BacterIal Genome-Scale Metabolic models for AppLied reverse ecoLogy
============

Python scripts and workflow by Matthew Jenior, University of Michigan, 2014 - 2016

---------------------------------------------------------------------------

The function of this package is to observe putative metabolic requirements, interaction, and activity of bacteria with the intention of inferring ecological interaction

KEGG reference files and reference creation script can be found in the support directory

Examples of input files for each program can be found in the examples directory

---------------------------------------------------------------------------

# Index:

**Section 1:**  Generating bipartite metabolic models used in conjunction with transcriptomics

**Section 2:**  Comparing bipartite metabolic model importance values between species

**Section 3:**  Appendix describing example files for each script

---------------------------------------------------------------------------

# Section 1 - bipartite_graph.py
Calculates relative importance of a given metabolite based on the expression of surrounding enzymes in a metabolic network

# Basic usage:
python bipartite_graph.py ko_expression.list

# Additional Options:
**Positional, required argument:**

expression_file

**Optional arguments:**

-h, --help		show this help message and exit

--name NAME		Organism or other name for KO+expression file (default is organism)

--iters ITERS		iterations for random distribution subsampling

---------------------------------------------------------------------------

##Under construction

# Section 2 - interact_bipartite.py

# Basic usage:
python interact_bipartite.py --files1 organism_1.scc.files --files2 organism_2.scc.files

# Additional Options:
**Required arguments:**

--files1 FILES1		 Directory of bipartite network output for first (meta)organism

--files2 FILES2		Directory of bipartite network output for second (meta)organism

**Optional arguments:**

-h, --help		show this help message and exit

--name1 NAME1		Name of first (meta)organism

--name2 NAME2		Name of second (meta)organism


---------------------------------------------------------------------------

# Section 3 - Appendix

**Sample files to be used as examples with each of the respective scripts**

ko_expression.tsv - list of KO codes and corresponding expression values for Clostridium difficile strain 630  
 
ko_expression.bipartite.files - Output of bipartite_graph.py for Clostridium difficile 630 during a mouse infection


