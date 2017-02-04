<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/mimosca_logo.png" title="MIMOSCA" alt="MIMOSCA" height=99 width=372>

# Multiple Input Multiple Output Single Cell Analysis

Named after <a href="https://en.wikibooks.org/wiki/Control_Systems/MIMO_Systems">MIMO control systems</a>, in which a combination of perturbations and measurements can help decode the dynamics of a system, this package was designed to assist those attempting to understand biological dynamics by designing, performing, and analyzing perturbation scRNA-seq experiments.

## Related Resources
* <a href="http://www.sciencedirect.com/science/article/pii/S0092867416316105">Our paper</a>
* <a href="https://groups.google.com/forum/#!forum/perturb-seq">Perturb-seq Google Forum</a>
* <a href="http://biorxiv.org/content/early/2016/12/12/093237">Chimera Correction</a>  and <a href="https://github.com/asncd/schimera">Code</a>
* <a href="http://www.clontech.com/US/Products/Genome_Editing/CRISPR_Cas9/Resources/Online_tools_for_guide_RNA_design">sgRNA Design Tools</a>
* Addgene Plasmids:  <a href="https://www.addgene.org/85967/">pBA439</a>
* <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90063">GEO Database link</a>

## Contents


* [FAQ](https://github.com/asncd/MIMOSCA/blob/master/README.md#faq)
* [Design of Experiments](https://github.com/asncd/MIMOSCA/blob/master/README.md#design-of-experiments--power-calculations)
* [Guide barcode, Cell barcode pairing](https://github.com/asncd/MIMOSCA/blob/master/README.md#guide-barcode-cell-barcode-pairing)
* [Computational Workflow](https://github.com/asncd/MIMOSCA/blob/master/README.md#computational-workflow)

Please let me know (here or in the Google Forum) if there are any areas that you'd like to see improved or new items to be added! 

## FAQ

* Q: How many perturbations can I do?
* A: It costs around $0.2/cell for commercial droplet scRNA-seq methods, it takes ~10 cells/perturbation to observe signature effects and ~100 cells/perturbation to see individual gene level effects robustly. Based on your budget, you can crunch the numbers (also see Design of Experiments section below)


## Design of Experiments & Power Calculations

<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/comp_knob.png" title="Experimental Design" alt="Experimental Design" height=360 width=432>

For a rough comparison of our pilot scRNA-seq to population RNA-seq of the same perturbation, see this <a href="https://github.com/asncd/MIMOSCA/blob/master/Power_Analysis_DOE/ost_ko_comparison.ipynb">iPython notebook</a>.

In designing Perturb-seq like experiments, there are a few key factors to keep in mind:

### Signatures vs. individual transcript-level phenotypes
Are you interested in broad transcriptional signatures or individual gene level differential expression? If the former, a rough approximation may be around 10 cells/perturbation. If the later, 100 or more cells may be required based on the effect size. 

A similar approximation for reads/cell would be a couple thousand for signatures and tens of thousands for gene-level.

### Library Size and Representation

As in any pooled screen, the representation of each perturbation in the library will vary. With genome wide CRISPR libraries the difference between the 10th and 90th percentile of a library is roughly 6-fold (Wang, 2013). Depending on how much a user wants to ensure every member of the library is represented, the cells/perturbation factor should be multiplied by an additional factor to reflect this variance.

### Using High MOI to infer genetic interactions

Our approach to use high MOI instead of either a single vector with multiple sgRNAs or vectors with different selection methods benefits from ease of implementation and the ability to represent a large diversity of combinations (only limited by the number of cells). 

However, challenges include a Poisson-like variance in the number of sgRNA/cell, sgRNA detection sensitivity, and the formation of PCR chimeras during the enrichment PCR procedure that can create misassignments. 

All three of these factors should be assessed in pilot experiments to troubleshoot. An example of such a pilot would look as follows (modified from the Drop-seq style species mixing experiments): 

<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/species_mix.png" title="Species Mix" alt="SMIX" height=431 width=365>

## Guide Barcode, Cell Barcode Pairing

<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/perturbseq_readdistribution.png" title="pseq_plasmid" alt="pseq_plasmid">

The distribution of reads going to the Perturb-seq vector (antiparallel) from 10X RNA-seq is shown above. Note that while the expression of the construct is comprable to that of a housekeeping gene, only a fraction of the reads overlap with the 18bp barcode (colored section in the coverage track). As such, it is advisable in most cases when you have a short barcode to perform enrichment PCR to obtain sensitive GBC/CBC pairing.

## Computational Workflow

<img src="https://github.com/asncd/MIMOSCA/blob/master/common_files/comp_flow2.png" title="Overview" alt="Overview" height=662 width=569>

### Inputs
* An expression matrix output by a high throughput scRNA-seq protocol (such as the <a href="http://mccarrolllab.com/dropseq/">Drop-seq</a> or <a href="https://support.10xgenomics.com/single-cell/software/pipelines/latest/what-is-cell-ranger">10X cellranger</a>)
* Guide barcode (GBC) PCR data to pair perturbations with cell barcodes (for certain applications this may be able to be directly obtained from the RNA-seq data
* A database of preassociated sgRNA/GBC pairs (either by Sanger sequencing or NGS)

### Intermediate Computation
* A simple fitness calculation is possible by determining the difference between the initial abundances of a GBC and how many cells it appeared in. 
* Guide barcodes and cell barcodes have to be paired accurately
* A Cell state classifier is defined on wildtype or control cells and then applied to all cells in an experiment. These classifications can used as outputs to be predicted (instead of gene expression) or as covariates in the model
* The linear model integrating all covariates (and interactions terms as desired) is fit. An EM-like approach filters cells that look much more like control cells than perturbed cells

### Outputs

* The regulatory coefficient obtained from the model are the most informative output giving an estimate of what extent each covariate (perturbation, cell state, pairwise interaction between perturbations, etc) impacted a given gene.
* Cell state effects are obtained by predicting the cell states based on the linear model instead of predicting gene expression
* Cell size effects (genes detected or transcripts detected) can be predicted as well


