BLASTats
========

Description
-----------

BLASTats BLASTs your protein sequence and tells you how it is distributed among <i>B. cereus</i> group.


Syntax
------

./bacilli.py [-v] [-l] [--identity] [--coverage] [--www protein_sequence]<br /><br />

-v: verbose mode<br /><br />
-l: lists precisely all hits in <i>B. cereus</i> group<br /><br />
--identity: manually set the identity threshold required between your query and a BLAST hit. Only hits above the threshold will be considered as potential functional homologues of your query protein.<br />Default value is 0.75.<br />For more information about functional homology infering, see:<br>Rost, B., Liu, J., Nair, R., Wrzeszczynski, K. O., & Ofran, Y. (2003). <i>Automatic prediction of protein function.</i> <b>Cellular and Molecular Life Sciences : CMLS</b>, 60(12), 2637–50.<br /><br />
--coverage: manually set the query coverage threshold required for a BLAST hit to be considered.<br /><br />
--www: followed by your protein sequence, this will BLAST it against NCBI's database<br />
