import os
import argparse
import modified_proteins_summ

parser = argparse.ArgumentParser(description='''Summarize post-translationally modified amino acid sites, peptides and proteins for Proteome Discoverer output''')

parser.add_argument('infile', metavar='-i', type=str, nargs='+', help='The input file should be PSMs table exported from Proteome Discoverer')

parser.add_argument('fasta', metavar='-f', type=str, nargs='+', help='Proteome database (FASTA) file used for the database search in Proteome Discoverer')

args = parser.parse_args()

def run_ptm_summ(infile, fasta):
    cmd = 'perl PTM_summarizer.pl ' + fasta + ' ' + infile
    os.system(cmd)
    print (cmd)

run_ptm_summ(args.infile[0], args.fasta[0])

