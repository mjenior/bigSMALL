bigsmall or *B*acter*I*al *G*enome-*S*cale *M*etabolic models for *A*pp*L*ied reverse eco*L*ogy
============

bigSMALL v1.0
Released: 12/01/2016

by
Matthew L. Jenior

Department of Microbiology & Immunology
University of Michigan
mljenior@umich.edu

When using, please cite:
Jenior, M.L., Leslie, J.L., Young, V.B., Schloss, P.D. (2016). Clostridium difficile colonizes alternative nutrient niches during infection across distinct murine gut environments. ...

Distributed under the GNU General Public License

---------------------------------------------------------------------------

Python scripts and workflow by Matthew Jenior, University of Michigan, 2014 - 2016

---------------------------------------------------------------------------

The function of this package is to infer putative metabolites most likely acquired for the environment based on transcriptomic data mapped to KEGG orthologs. Bigsmall generates a bipartite metabolic based on the reaction data associated with each KEGG ortholog and then integrates transcript abundances to predict demand for metabolites based on the transciption of adjacent enzyme nodes. Monte Carlo simulation is also applied to create a standard of comparison that reflects random noise, where highly dinstinct measurments are more likely to be not by chance.

KEGG reference files and reference creation script can be found in the support directory

Examples of input files for each program can be found in the examples directory

---------------------------------------------------------------------------

# bigsmall.py
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

* A good starting point is to use x10 the number of genes for --iters

---------------------------------------------------------------------------

# Supporting files

**Sample files to be used as examples with each of the respective scripts**

ko_expression.tsv - list of KO codes and corresponding expression values for Clostridium difficile strain 630  
 
ko_expression.bipartite.files - Output of bipartite_graph.py for Clostridium difficile 630 during a mouse infection

---------------------------------------------------------------------------

# Citations

McGill, R., Tukey, J. W. & Larsen, W. A. (1978). Variations of Box Plots. The American Statistician 32, 12–16.

Ogata, H. et al. (1999). KEGG: Kyoto encyclopedia of genes and genomes. 27, 29–34.

Patil, K. R. & Nielsen, J. (2004). Uncovering transcriptional regulation of metabolism by using metabolic network topology. PNAS 102 (8), 2685–2689.
