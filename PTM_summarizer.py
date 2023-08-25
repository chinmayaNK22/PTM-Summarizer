from itertools import islice
from read_fasta_file_v2 import readfasta
import os
import modified_proteins_summ
import argparse

parser = argparse.ArgumentParser(description='''Summarize post-translationally modified amino acid sites, peptides and proteins for Proteome Discoverer output''')

parser.add_argument('infile', metavar='-i', type=str, nargs='+', help='The input file should be PSMs table exported from Proteome Discoverer')

parser.add_argument('fasta', metavar='-f', type=str, nargs='+', help='Proteome database (FASTA) file used for the database search in Proteome Discoverer')

args = parser.parse_args()

def get_header_idx(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            if '"' in split_i[0]:
                try:
                    pep = split_i.index('"Annotated Sequence"')
                except:
                    pep = split_i.index('"Sequence"')
                pro = split_i.index('"Master Protein Accessions"')
                mz = split_i.index('"m/z [Da]"')
                scan = split_i.index('"First Scan"')
                try:
                    prob_score = split_i.index('"ptmRS: Best Site Probabilities"')
                except :
                    raise ("ERORR: There is no ptmRS: Best Site Probabilities column present in the file")
            else:
                try:
                    pep = split_i.index("Annotated Sequence")
                except:
                    pep = split_i.index('Sequence')
                pro = split_i.index('Master Protein Accessions')
                mz = split_i.index('m/z [Da]')
                scan = split_i.index('First Scan')
                try:
                    prob_score = split_i.index('ptmRS: Best Site Probabilities')
                except :
                    raise ("ERORR: There is no ptmRS: Best Site Probabilities column present in the file")                

            return pep, pro, mz, scan, prob_score

def parse_ptmRS_score(instring):
    #print (instring)
    try:
        aa = instring[0]
        pos = instring.split('(')[0].lstrip(aa)
        mod = instring.split('(')[1].split(')')[0]
        score = instring.split(':')[1].strip()
        return aa, pos, mod, score
    except:
        print (instring)

def parse_psm_file(infile):
    a = get_header_idx(infile)
    mod_scores = {}
    modifications = {}
    with open(infile) as file:
        for i in islice(file,1,None):
            split_i = i.rstrip().split('\t')
            pep = split_i[a[0]].strip('"').split('.')[1]
            pro = split_i[a[1]].strip('"')
            mz = split_i[a[2]].strip('"')
            scan = split_i[a[3]].strip('"')
            if a[4] <= len(split_i):
                ptmrs_score = split_i[a[4]].strip('"')
                if len(ptmrs_score) != 0:
                    new_mod_score = {}
                    if ';' in ptmrs_score:
                        for mod_score in ptmrs_score.split(';'):
                            AA, POS, MOD, SCORE = parse_ptmRS_score(mod_score.strip())
                            if float(SCORE) > 75.0:
                                new_mod_score[AA + POS + '(' + MOD + '): ' + SCORE] = mod_score
                                if MOD + '_' + AA not in modifications:
                                    modifications[MOD + '_' + AA] = [split_i]
                                else:
                                    modifications[MOD + '_' + AA].append(split_i)
                                #mod_scores[pep + '@' + pro + '@' + mz + '@' + scan + '@' + AA + '@' + POS + '@' + MOD + '@' + SCORE] = [split_i]
           
                    elif ptmrs_score.strip() != 'Too many isoforms':
                        AA, POS, MOD, SCORE = parse_ptmRS_score(ptmrs_score.strip())
                        if float(SCORE) > 75.0:
                            new_mod_score[AA + POS + '(' + MOD + '): ' + SCORE] = ptmrs_score
                            if MOD + '_' + AA not in modifications:
                                modifications[MOD + '_' + AA] = [split_i]
                            else:
                                modifications[MOD + '_' + AA].append(split_i)
                            #mod_scores[pep + '@' + pro + '@' + mz + '@' + scan + '@' + AA + '@' + POS + '@' + MOD + '@' + SCORE] = [split_i]

                    final_score = {k:v for k, v in new_mod_score.items() if len(k) != 0}
                    score_col = ';'.join(list(final_score))
                    if len(score_col) != 0:
                        mod_scores[pep + '@' + pro + '@' + mz + '@' + score_col] = [split_i]
                        
    return mod_scores, modifications

def parse_acc(header):
    splitters = [' ','|']
    acc_sep = []
    for split in splitters:
        acc_sep.append(header.index(split))
        
    return sorted(acc_sep)[0]
    
def map_to_protein(indict, infasta):
    output = {}
    for rows in readfasta(infasta).read():
        header = rows[0]
        seq = rows[1]
        splitter = parse_acc(header)
        acc = header[0:splitter]
        for keys, values in indict.items():
            mod_peps = keys.split('@')
            if ';' in mod_peps[1]:
                for pro in mod_peps[1].split(';'):
                    if acc == pro:
                        if mod_peps[0].upper() in seq:
                            pep_pos = []
                            pro_pos = []
                            aa = []
                            mods = []
                            for best_ptmrs in mod_peps[3].split(';'):
                                AA, POS, MOD, SCORE = parse_ptmRS_score(mod_peps[3])
                                if seq[seq.index(mod_peps[0].upper())+ (int(POS)-1)] == AA:
                                    pep_pos.append(str(POS))
                                    pro_pos.append(str(seq.index(mod_peps[0].upper())+ int(POS)))
                                    aa.append(seq[seq.index(mod_peps[0].upper())+ (int(POS)-1)])
                                    mods.append(MOD)
                                    
                            output[mod_peps[0] + '@' + mod_peps[1] + '@' +  ';'.join(mods) + '@' + ';'.join(aa) + '@' + ';'.join(pep_pos) + '@' + ';'.join(pro_pos)] = [keys]
                            
                            #if seq[seq.index(mod_peps[0].upper())+ (int(mod_peps[5])-1)] == AA:
                            #    output[mod_peps[0] +'@'+ mod_peps[1] +'@'+ mod_peps[-2]+'@'+ mod_peps[4]+'@'+ (mod_peps[-2] + '_' + mod_peps[4]) +'@'+ str(seq.index(mod_peps[0].upper()) + (int(mod_peps[5])-1))] = [keys]
                            #else:
                            #    print ("Modified amino acid sequence ", mod_peps[0], " is not present in protein ", acc)
            else:
                if acc == mod_peps[1]:
                    if mod_peps[0].upper() in seq:
                        pep_pos = []
                        pro_pos = []
                        aa = []
                        mods = []
                        for best_ptmrs in mod_peps[3].split(';'):
                            AA, POS, MOD, SCORE = parse_ptmRS_score(best_ptmrs)
                            if seq[seq.index(mod_peps[0].upper())+ (int(POS)-1)] == AA:
                                pep_pos.append(str(POS))
                                pro_pos.append(str(seq.index(mod_peps[0].upper())+ int(POS)))
                                aa.append(seq[seq.index(mod_peps[0].upper())+ (int(POS)-1)])
                                mods.append(MOD)
                                #print (mod_peps[0], mod_peps[-1], POS, seq.index(mod_peps[0].upper())+ (int(POS)-1), seq[seq.index(mod_peps[0].upper())+ (int(POS)-1)])
                                          
                        output[mod_peps[0] + '@' + mod_peps[1] + '@' +  ';'.join(mods) + '@' + ';'.join(aa) + '@' + ';'.join(pep_pos) + '@' + ';'.join(pro_pos)] = [keys]
                

    outputs = []
    for k, v in output.items():
        outputs.append(k.split('@'))

    return outputs

def write_to_file(outlist, infile):
    outfile = "{0}_UniqueProteinSite.txt".format(infile.rstrip('txt').rstrip('.'))
    with open(outfile, 'w') as outf:
        outf.write('Peptide\tProtein\tModification\tAmino_Acid\tModified_Site (Peptide)\tPTM_Site(Protein)\n')
        outf.writelines('\t'.join(i) + '\n' for i in outlist)

def write_mod_files(mod_dicts, infile):
    header = open(infile).readline().rstrip().split('\t')
    for k, v in mod_dicts.items():
        outfile = infile.rstrip('txt').rstrip('.') + '_' + k + '.txt'
        with open(outfile, 'w') as outf:
            outf.write('\t'.join(header) + '\n')
            outf.writelines('\t'.join(i) + '\n' for i in v)

def summarize_ptm(infile, infasta):
    modified_psms, mod_psms = parse_psm_file(os.path.join(infile))

    ### Write modification specific PSMs to seperate files based on ptmRS probability threshold
    write_mod_files(mod_psms, infile)

    ### Map modified peptides to proteins and fetch modification site at protein level
    output = map_to_protein(modified_psms, os.path.join(infasta))
    
    write_to_file(output, infile)
    print (len(output))
    if len(output) != 0:
        modified_proteins_summ.summarize_ptms(infile, output)
    else:
        print ("No modifications were found with ptmRs probability score > 75")


summarize_ptm(args.infile[0], args.fasta[0])
