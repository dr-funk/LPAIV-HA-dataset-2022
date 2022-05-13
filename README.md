# LPAIV-HA-dataset-2022
Scripts to analyse the dataset from Funk et al. 2022 

These scripts are meant to be used with a sorted dataset like the one provided alongside the article. Data has to be sorted by subtype as well as splitting variable(s) (region and species in this case). The script will produce graphs taking only the terrestrial poultry and wild aquatic bird files into account.

For all scripts the folder containing the data, i.e., the subfolders split by region and species, has to be provided. The scripts require the data to be organized in the same way as the dataset from the article, i.e.,
```
<data folder>
   └── by_region-species
         ├── america-ans_cha
         │     └── fasta files in 'america-ans_cha_hX_unique_CSregion.fasta' format 
         ├── america-poultry
         │     └── fasta files in 'america-poultry_hX_unique_CSregion.fasta' format 
         ├── eurasia-ans_cha
         │     └── fasta files in 'eurasia-ans_cha_hX_unique_CSregion.fasta' format 
         └── eurasia-poultry
               └── fasta files in 'eurasia-poultry_hX_unique_CSregion.fasta' format 
 ```
