from itertools import islice
from collections import Counter
##import argparse
##
##parser = argparse.ArgumentParser(description='''Summarize the modification found in all the modified proteins''')
##
##parser.add_argument('infile', metavar='-i', type=str, nargs='+', help='Output file from PTM-Summarizer.pl')
##
##args = parser.parse_args()

def modification_summ(inlist):
    dicts = {}
    dicts1 = {}
    dicts2 = {}
    summary = {}
    for i in inlist:
        print (i)
        if i[1] not in dicts:
            dicts[i[1]] = [i[4]]
            #dicts1[split_i[1]] = [split_i[2] + '\t' + split_i[3] + '\t' + split_i[4] + '\t' + split_i[5].rstrip()]
        else:
            dicts[i[1]].append(i[4])
            #dicts1[split_i[1]].append(split_i[2] + '\t' + split_i[3] + '\t' + split_i[4] + '\t' + split_i[5].rstrip())
        if i[-2] not in dicts2:
            summary[i[-2]] = [i[1]] 
            dicts2[i[-2]] = [1]
        else:
            summary[i[-2]].append(i[1])
            dicts2[i[-2]].append(1)
            
    return summary, dicts, dicts2

##if args.infile[0].split('_')[-1] == '_UniqueProteinSite.txt':
##    modification_summ(args.infile[0])
##else:
##    print ('ERROR: The input file in not correct.')
###file = 'SARS-CoV-2_Multi-PTM_Unique_Modification_Sites_062821.txt'

def summarize_ptms(infile, inlist):
    summary, proteins, modifications = modification_summ(inlist)
    summ_output = [[k.split('_')[0], k.split('_')[1], str(len(v))] for k, v in summary.items()]
    summ_outfile = '{0}_Summary.txt'.format(infile.rstrip('txt').rstrip('.'))
    with open(summ_outfile, 'w') as sumf:
        sumf.write('Modification\tAmino Acid\tTotal No. of modified sites\n')
        sumf.writelines('\t'.join(i) + '\n' for i in summ_output)

    output = []
    for k, v in proteins.items():
        dicts_vs = {}
        for x in v:
            if x not in dicts_vs:
                dicts_vs[x] = [1]
            else:
                dicts_vs[x].append(1)

        output_dict = {}
        for w, z in modifications.items():
            if w in dicts_vs:
                if k not in output_dict:
                    output_dict[k] = [str(len(dicts_vs[w]))]
                else:
                    output_dict[k].append(str(len(dicts_vs[w])))
            else:
                if k not in output_dict:
                    output_dict[k] = [str(0)]
                else:
                    output_dict[k].append(str(0))
        for l, m in output_dict.items():
            output.append([l] + m + [str(sum([int(v) for v in m]))])

    header = [['Accession'] + list(modifications) + ['Total No. of Modified Sites']]
    #print (header)
    outfile = '{0}_Modified_Proteins.txt'.format(infile.rstrip('txt').rstrip('.'))
    with open(outfile, 'w') as outf:
        outf.writelines('\t'.join(i) + '\n' for i in header)
        outf.writelines('\t'.join(i) + '\n' for i in output)
