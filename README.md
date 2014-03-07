BLASTats
========

Description
-----------

BLASTats BLASTs your protein sequence and tells you how it is distributed among <i>B. cereus</i> group.


Syntax
------

./bacilli.py [-v] [-l] [--www protein_sequence]

-v: verbose mode<br />
-l: lists precisely all hits in <i>B. cereus</i> group<br />
--www: followed by your protein sequence, this will BLAST it against NCBI's database<br />
--identity: manually set the identity threshold required between your query and a BLAST hit. Only hits above the threshold will be considered as potential functional homologues of your query protein. Default value is 0.75. For more information on functional homology infering, see: Rost, B., Liu, J., Nair, R., Wrzeszczynski, K. O., & Ofran, Y. (2003). <i>Automatic prediction of protein function.</i> <b>Cellular and Molecular Life Sciences : CMLS</b>, 60(12), 2637–50.<br />
--coverage: manually set the query coverage threshold required for a BLAST hit to be considered.<br />