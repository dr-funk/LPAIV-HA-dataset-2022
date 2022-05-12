# all_analyses_interactive.py
This script generates an interactive figure like Supplementary Figure 3 from the article. It automatically loops through all analyses and ways to split the data and outputs a single HTML file with interactive tabs to look through the data. The following analyses are performed by default, but the set-up is modular allowing to add or remove analyses:
  * Counting the number of adenines in each sequence
  * Finding the length of the longest adenine stretch in each sequence
  * Counting the number of purines in each sequence
  * Finding the length of the longest purine stretch in each sequence
  * Counting the minimal number of substitutions required to reach a (non-consecutive) tribasic cleavage site for each sequence
  * Counting the minimal number of substitutions required to have an adenine stretch of at least a threshold (varied from 3-10) for each sequence

The script requires [Biopython](https://anaconda.org/conda-forge/biopython), [pandas](https://anaconda.org/conda-forge/pandas), [Plotly](https://anaconda.org/plotly/plotly) and [jinja2](https://anaconda.org/anaconda/jinja2).

The folder containing the sorted data needs to be provided, there are two ways to do this, either on the command line using
`all_analyses_interactive.py <path to data folder>`
or by writing the path to the data folder on line 21 in which case no command line arguments are required.

The jinja template needs to be placed into the data folder or its enclosing folder needs to be specified on line 24.
