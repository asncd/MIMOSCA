This folder contains iPython notebooks that process Drop-seq output (examples included in the data subfolder) with **scPerseus** (single cell Perseus) and **Theseus** and a python script that processes 10X output run using cellranger 1.2+ with **Theseus**. The names [Bellerophon](http://comp-bio.anu.edu.au/Bellerophon/doc/doc.html) and [Perseus](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-38) were taken already in the context of PCR chimeras.

![alt tag](http://www.greekmythology.com/images/mythology/theseus_adventures_78.jpg)

**Theseus** is a relatively brute force approach to tackling the chimera problem by filtering molecules that, for a given cell barcode and UMI pair, have low read abundance.

<img src="http://i.imgur.com/olgUb2b.jpg" alt="Perseus" width="307" height="203">

**scPerseus**, is a slighlty more sophisticated option, in which the user specifies labels for cells and a set of genes that should be relatively unique for that population. A random forest classifier is trained on the dataset using a diverse set of features. A filter is applied using the out of bag probability estimates to remove chimeric molecules. 


Usage for 10x **Theseus** script:


python theseus_10x.py  FULLPATH2H5.h5 OUTPUTPATH TPT_THRESHOLD NUMBER_OF_CELLS
