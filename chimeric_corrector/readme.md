This folder contains iPython notebooks that process dropseq output (examples included in the data subfolder) with Perseus and Theseus (the name [Bellerophon]: http://comp-bio.anu.edu.au/Bellerophon/doc/doc.html was taken already)


![alt tag](http://www.greekmythology.com/images/mythology/theseus_adventures_78.jpg)

Theseus is a relatively brute force approach to tackling the chimera problem by filtering molecules that, for a given cell barcode and UMI pair, have low read abundance.






![alt tag2](http://i.imgur.com/olgUb2b.jpg)

Perseus, is a slighlty more sophisticated option, in which the user specifies labels for cells and a set of genes that should be relatively unique for that populaiton. A random forest classifier is trained on the dataset using a diverse set of features. A filter is applied using the out of bag probability estimates to remove chimeric molecules. 
