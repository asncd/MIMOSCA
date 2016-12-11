# Theseus
<img align="left" src="http://www.greekmythology.com/images/mythology/theseus_adventures_78.jpg" title="Theseus defeating the Minatour" alt="Theseus defeating the Minatour" width="250" height="247">

**Theseus** is a relatively brute force approach to tackling the PCR chimera problem in scRNA-seq by filtering molecules that, for a given cell barcode and UMI pair, have low read abundance. We call this TPT normalization and is defined as follows. For a given pair of cell barcode and UMI, *u*, and a particular transcript *i*:

TPT<sub>i,u</sub> = (reads<sub>i,u</sub>)/sum(reads<sub>j,u</sub>)

Where the denominator is summed over all transcripts *j* sharing the same *u*.

The user can define a threshold to filter molecules with low TPT and obtain a new expression matrix.


# Yorimasa

<img align="right" src="https://data.ukiyo-e.org/famsf/images/6340304231510089.jpg" title="Yorimasa defeating the Nue" alt="Yorimasa defeating the Nue" width="307" height="203">

**Yorimasa**, is a slighlty more sophisticated option, in which the user specifies labels for cells and a set of genes that should be relatively unique for that population. A random forest classifier is trained on the dataset using a diverse set of features. A filter is applied using the out-of-bag probability estimates to remove chimeric molecules. 

The features used may include:
* Reads per molecule
* TPT
* Relative abundance of the cell barcode
* Relative abundance of the transcript
* Transcript GC content
* Transcript length

# Folder Contents

This folder contains iPython notebooks that process Drop-seq output (examples are included in the data subfolder along with precomputed gene level features) with **Yorimasa** and **Theseus** and a python script that processes 10X output run using cellranger 1.2+ with **Theseus**.  The names [Bellerophon](http://comp-bio.anu.edu.au/Bellerophon/doc/doc.html) and [Perseus](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-38) were taken already in the context of PCR chimeras.

Usage for 10x **Theseus** script:

```python
python theseus_10x.py  FULLPATH2H5.h5 OUTPUTPATH TPT_THRESHOLD NUMBER_OF_CELLS
```
