# PTM-Summarizer
The modified peptides and sites identified from proteome database search using Proteome Discoverer can be summarized using PTM-Summarizer. This tool requires Peptide Spectrum Match (PSM) output file from PTM search with "ptmRS: Best Site Probabilities" column and the proteome database (FASTA format) used for the search.


**How to run PTM-Summarizer:**
```
PTM-Summarizer>python HumanRefSeq_GCF_000001405.39_GRCh38.p13_protein_refseq109_formated.fasta SARS_CoV2_Mulit-PTM_PSMs_PTMs.txt
```
The output of this tool consists of;
1. List of modified site information for each modification in each protein (Infile_UniqueProteinSite.txt)
2. Summary of the modification specific sites (Infile_summary.txt)
3. Input file for Upset plot (ForUpset.txt)
4. Modification and site specific test file (ex: Acetyl_K.txt, Crotonyl_K.txt etc)

For any enquiries please contact: beheras40@gmail.com and sandeep.kolya@gmail.com or chinnu.kemmaai@gmail.com
