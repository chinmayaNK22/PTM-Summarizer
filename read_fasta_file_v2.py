from itertools import groupby


class readfasta(object):

    def __init__(self, infile):
        self.infile = infile

    def read(self):
        file_format = {'fasta':1, 'faa':2, 'fna':3}
        if self.infile.split('.')[-1].lower() in file_format:
            read_file = open(self.infile)
            faiter = (x[1] for x in groupby(read_file, lambda line: line[0] == ">"))
            for header in faiter:
                header = next(header)[1:].strip()
                seq = "".join(s.strip() for s in next(faiter))
                yield header, seq

        else:
            filetype = self.infile.split('.')[-1]
            raise Exception(f'Can only read FASTA file format with file extensions fasta, faa or fna')


class decoy_pro(object):
    def __init__(self, protein):
        self.protein = protein

    def Reverse(self):
        pro = [self.protein[p] for p in range(len(self.protein)) if p+1 == len(self.protein)] + [self.protein[p] for p in range(len(self.protein)) if p+1 != len(self.protein)]
        pro.reverse()
        return ''.join(pro)
