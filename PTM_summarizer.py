import os
import argparse

parser = argparse.ArgumentParser(description='''SSPE: Species Specific Peptide Extractor. Peptides list from the species specific input fasta file (Proteome) will be
                                            generated and only those peptide sequences specific to a input fasta file (i.e. Species proteome) will be considered when compared with the other input proteome databases''')

parser.add_argument('infile', metavar='-i', type=str, nargs='+', help='The input file should be PSMs table exported from Proteome Discoverer')

parser.add_argument('fasta', metavar='-fa', type=str, nargs='+', help='Proteome database (FASTA) file used for the database search in Proteome Discoverer')

args = parser.parse_args()

cmd = 'perl PTM_summarizer.pl ' + args.fasta[0] + ' ' + args.infile[0]
os.system(cmd)
print (cmd)
