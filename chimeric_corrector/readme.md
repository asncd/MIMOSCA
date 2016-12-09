This folder contains iPython notebooks that process Drop-seq output (examples included in the data subfolder) with **Perseus** and **Theseus** and a python script that processes 10X output run using cellranger 1.2+ with **Theseus**. The name [Bellerophon](http://comp-bio.anu.edu.au/Bellerophon/doc/doc.html) was taken already.

![alt tag](http://www.greekmythology.com/images/mythology/theseus_adventures_78.jpg)

**Theseus** is a relatively brute force approach to tackling the chimera problem by filtering molecules that, for a given cell barcode and UMI pair, have low read abundance.

<img src="http://i.imgur.com/olgUb2b.jpg" alt="Perseus" width="307" height="203">

**Perseus**, is a slighlty more sophisticated option, in which the user specifies labels for cells and a set of genes that should be relatively unique for that population. A random forest classifier is trained on the dataset using a diverse set of features. A filter is applied using the out of bag probability estimates to remove chimeric molecules. 
