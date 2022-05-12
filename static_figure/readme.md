# Scripts to generate static pdf files

Both scripts provided here can be used to generate static pdf graphs similar to those in the article. The scripts requires [Biopython](https://anaconda.org/conda-forge/biopython), [pandas](https://anaconda.org/conda-forge/pandas), [Matplotlib](https://anaconda.org/conda-forge/matplotlib) and [Seaborn](https://anaconda.org/anaconda/seaborn).

## muts_to_basic.py

This script finds the minimum number of substitutions required to reach a (non-)consecutive tri- or tetrabasic site for each sequence of a dataset. There are several parameters that can be changed:
  * Line 29: if set to True, the basic amino acids need to be consecutive else. This only impacts tribasic cleavage sites
  * Line 30: if set to True, a tribasic site is required, if set to False a tetrabasic site is required
  * Line 31: the target motif

The folder containing the sorted data needs to be provided as well as the second variable by which to split the data (region, species, region-species or none). This can be done in two ways, either via command line arguments `muts_to_basic.py <path to data folder> <variable2>` of by writing the path to the data folder on line 15 and the variable on line 21. By default line 21 is set to none, so that if only the path to the data folder is provided as command line argument, the analysis is run on all sequences. If provided command line arguments will override the hard-coded arguments.

## muts_to_astretch.py

This script finds the minimum number of substitutions required to reach an A stretch of over a certain threshold (varied from 3-10) for each sequence of a dataset. 

The folder containing the sorted data needs to be provided as well as the second variable by which to split the data (region, species, region-species or none). This can be done in two ways, either via command line arguments `muts_to_astretch.py <path to data folder> <variable2>` of by writing the path to the data folder on line 15 and the variable on line 21. By default line 21 is set to none, so that if only the path to the data folder is provided as command line argument, the analysis is run on all sequences. If provided command line arguments will override the hard-coded arguments.
