from itertools import islice
from collections import Counter

def for_upsetplot(inlist, infile):
    ''' Generate a input for UpSet plot on number of modified protein sites '''
    mods = {}
    dicts = {}
    for i in inlist:
        for idx, j in enumerate(i[-1].split(';')):
            pep_pos = str(i[-1].split(';')[idx])
            aa = i[-3].split(';')[idx]
            mod = i[2].split(';')[idx]
            if i[1] + '_' + pep_pos not in dicts:
                dicts[i[1] + '_' + pep_pos] = [mod + '_' + aa]
            else:
                dicts[i[1] + '_' + pep_pos].append(mod + '_' + aa)

            mods[mod + '_' + aa] = mod

    outfile = "{0}_UpSet.txt".format(infile.rstrip('txt').rstrip('.'))
    header = ['Protein_Site'] + list(mods)
    with open(outfile, 'w') as outf:
        outf.write('\t'.join(header) + '\n')
        for k, v in dicts.items():
            nr_mod = {j:j for j in v}
            val = []
            for mod in list(mods):
                if mod in nr_mod:
                    val.append(str(1))
                else:
                    val.append(str(0))

            
            outf.write(k + '\t' + '\t'.join(val) + '\n')
        
def modification_summ(inlist):
    ''' Summarize the information on number of modifications and protein sites '''
    dicts = {}
    dicts1 = {}
    dicts2 = {}
    summary = {}
    for i in inlist:
        for idx, j in enumerate(i[-1].split(';')):
            aa = i[-3].split(';')
            mods = i[2].split(';')
            mod_aa = mods[idx] + '_' + aa[idx]
            
        if i[1] not in dicts:
            dicts[i[1]] = [mod_aa]  
        else:
            dicts[i[1]].append(mod_aa)
            
        if mod_aa not in dicts2:
            summary[mod_aa] = [i[1]] 
            dicts2[mod_aa] = [1]
        else:
            summary[mod_aa].append(i[1])
            dicts2[mod_aa].append(1)
            
    return summary, dicts, dicts2


def summarize_ptms(infile, inlist):
    
    summary, proteins, modifications = modification_summ(inlist)
    
    for_upsetplot(inlist, infile)
    
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
