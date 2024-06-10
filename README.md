# eMMA
Functional Families generation using embedding distance matrices

Overview
---------

New CATH-Gemma algorithm to generate trees of relationships between sequences and functions.

The key change allows to use distances from pLM embeddings or structural distances instead of HMM-vs-HMM comparisons. 

### Main features

- Revised protocol to use MMseqs2 instead of CD-HIT. 
- Python CLI generating SGE or local jobs
- embedding distances or 1/bitscore distances from Foldseek as data source for functional relationships
- Faster, low memory footprint. (i.e. For the HUPS Superfamily (3.40.50.620) 22 hours to 6 hours).

This repo is part of the FunFams pipeline as an intermediate step before FunFHMMER. 

The eMMA version of FunFHMMER can be found at [funfhmmer-emma]([https://github.com/UCL/cath-funfhmmer/tree/funfhmmer-emma](https://github.com/UCLOrengoGroup/eMMA-FunFHMMER))

See the GeMMA [Wiki](https://github.com/UCL/cath-gemma/wiki) for documentation on GEMMA and check out the step-by-step walkthrough [here](https://github.com/UCLOrengoGroup/eMMA/blob/main/step-by-step-walkthrough.md)
