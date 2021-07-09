# PTM-Summarizer
The modified peptides and sites identified from proteome database search using Proteome Discoverer can be summarized using PTM-Summarizer. This tool requires Peptide Spectrum Match (PSM) output file from PTM search with "ptmRS: Best Site Probabilities" column and the proteome database (FASTA format) used for the search.


**How to run PTM-Summarizer:**
```
PTM-Summarizer>perl PTM_ProAdvanced.pl
```
The output of this tool consists of;
1. List of modified site information for each modification in each protein (UniqueProteinSite.txt)
2. Summary of the modification specific sites (summary.txt)
3. Input file for Upset plot (ForUpset3.txt)
4. Modification and site specific test file (ex: Acetyl_K.txt, Crotonyl_K.txt etc)
