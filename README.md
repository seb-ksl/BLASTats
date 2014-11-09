BLASTats
========

Description
-----------

BLASTats is a pipeline written in Python, using the BioPython module. It allows getting a one-click idea of a protein distribution among species of interest, by BLASTing a given protein sequence, parsing the results and building a phylogenetic tree of homologues.

Find out more here: http://gelis.ch/programs/blastats/

Changelog v 2.0
---------------

- Lengthened pipeline: adding homologues alignment with ClustalOmega and building phylogenetic tree with FastTree
- Fixed formula for query coverage calculation
- Refined regex to solve bug when only 1 genome available for a specie
- Enhanced output