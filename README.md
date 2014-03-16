BLASTats
========

Description
-----------

BLASTats helps you see how your protein of interest is distributed among <i>B. cereus</i> group.
This program must be run in a console (Linux), terminal (Mac), or cmd (Windows).

Syntax
------

./blastats.py(c) [-v] [-l] [--identity] [--coverage] [--www protein_sequence]

Parameters
----------

* -v: verbose mode
* -l: lists precisely all hits in <i>B. cereus</i> group
* --identity: manually set the identity threshold required between your query and a BLAST hit. Only hits above the threshold will be considered as potential functional homologues of your query protein.<br />Default value is 0.75.<br />For more information about functional homology infering, see:<br>Rost, B., Liu, J., Nair, R., Wrzeszczynski, K. O., & Ofran, Y. (2003). <i>Automatic prediction of protein function.</i> <b>Cellular and Molecular Life Sciences : CMLS</b>, 60(12), 2637–50.
* --coverage: manually set the query coverage threshold required for a BLAST hit to be considered.
* --www: followed by your protein sequence, this will BLAST it against NCBI's database.
