#!/usr/bin/python3
# my functions - 10/20/2025
# Hee Jin Cho

def myfunc1():
    '''this is a docstring'''
    # this is a comment
    pass

def helloworld():
    print('hello, world')

def echo(whatever):
    print(whatever)

def addsub_x(x, y): # not print due to the return
    a = x + y
    b = x - y
    return (a, b)
    print('%s + %s = %s' % (x, y, a))
    print('%s - %s = %s' % (x, y, b))


def addsub(x, y):
    a = x + y
    b = x - y
    print('%s + %s = %s' % (x, y, a))
    print('%s - %s = %s' % (x, y, b))
    return (a, b)

def welcome(x, y = 'hello'):
    print(f'{y}, {x}')

def welcome_bst(x = 'bst', y = 'hello'):
    print(f'{y}, {x}')

#wrong example
#def welcome_bst2(x = 'bst', y):
#    print(f'{y}, {x}')

def echo_list(*args):
    print(args)
    print([x for x in args])

def lang_ver(lang, *args):
    ver_str = ', '.join([str(ver) for ver in args])
    res_str = '%s versions: %s' % (lang, ver_str)
    return res_str

def reverse_complement(orig_seq):
    nt_pairH = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    temp_seq = orig_seq[::-1]
    rev_seq = ''.join([nt_pairH.get(nt, 'X') for nt in temp_seq])
    return rev_seq

def reverse_complement_list(*args):
    nt_pairH = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    res_list = []
    for orig_seq in args:
        temp_seq = orig_seq[::-1]
        rev_seq = ''.join([nt_pairH.get(nt, 'X') for nt in temp_seq])
        res_list.append(rev_seq)
    return res_list

def dict_return(**kwargs):
    return kwargs

def add2num(x,y):
    z = float(x) + float(y)
    return z

def add2num_rev(x,y):
    try:
        z = float(x) + float(y)
        return z
    except:
        print('only int/float are allowed')
    finally:
        print('this function is working')

def oddeven(x):
    if not isinstance(x, (int)):
        raise Exception(f'Integers are only accepted. Your input is {type(x)}')
    elif x % 2 == 1:
        print('x is odd')
    else:
        print('x is even')

class MySeq:
    '''this is from the textbook'''
    def __init__(self, seq, seq_type = "DNA"):
        self.seq = seq
        self.seq_type = seq_type
    def print_sequence(self):
        print('Sequence: ' + self.seq)
    def get_seq_biotype(self):
        return self.seq_type
    def show_info_seq(self):
        print('Sequnece: ' + self.seq + " biotype: " + self.seq_type)
    def count_occurrences(self, seq_search):
        return self.seq.count(seq_search)

class rcClass:
    nt_pairH = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    def __init__(self, orig_seq):
        self.orig_seq = orig_seq
        self.temp_seq = self.orig_seq[::-1]
        self.rev_seq = ''.join([rcClass.nt_pairH.get(nt, 'X') for nt in self.temp_seq])
    def reverse_complement(self):
        return self.rev_seq
    def nt_count(self, count_nt):
        return self.orig_seq.count(count_nt)
    def rc_nt_count(self, count_nt):
        return self.rev_seq.count(count_nt)

class rcClassChild(rcClass):
    def nt_count(self, count_nt): #overriding
        if self.orig_seq.count(count_nt) == 0:
            print('"%s" is not in "%s"' % (count_nt, self.orig_seq))
            return 0
        else:
            return self.orig_seq.count(count_nt)
    def nt_ratio1(self, count_nt):
        return self.nt_count(count_nt)/len(self.orig_seq)
    def nt_ratio2(self, count_nt):
        return super().nt_count(count_nt)/len(self.orig_seq)
    def rc_nt_ratio(self, count_nt):
        return super().rc_nt_count(count_nt)/len(self.rev_seq)
#    def nt_ratio(self, nt):
#        return self.nt_count(nt)/len(self.rev_seq)

class rcClassChild2(rcClassChild):
    def __init__(self, orig_seq):
        self.orig_seq = orig_seq.upper()
        self.temp_seq = self.orig_seq[::-1]
        self.rev_seq = ''.join([rcClass.nt_pairH.get(nt, 'X') for nt in self.temp_seq])

def fileCapitalize(inFilename, outFilename):
    in_file = open(inFilename,'r')
    out_file = open(outFilename, 'w')
    for line in in_file:
        out_file.write(line[:-1].upper() + '\n')
        #out_file.write(line.upper())
    in_file.close()
    out_file.close()

def translate_codon(cod):
    """ Translates a codon into an amino acid using an internal dictionary with the standard genetic code."""
    tc = {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TGT": "C", "TGC": "C", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TTT": "F", "TTC": "F", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "CAT": "H", "CAC": "H", "ATA": "I", "ATT": "I", "ATC": "I", "AAA": "K", "AAG": "K", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "ATG": "M", "AAT": "N", "AAC": "N", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TGG": "W", "TAT": "Y", "TAC": "Y", "TAA": "_", "TAG": "_", "TGA": "_"}
    if cod in tc:
        return tc[cod]
    else:
        return None

def reading_frames(dna_seq):
    """ Computes the six reading frames of a DNA sequence including the reverse complement."""
    assert valid_dna(dna_seq), 'Your sequence is not valid'
    res_orig = [translate_seq(dna_seq, i) for i in range(0,3)]
    res_rc = [translate_seq(reverse_complement(dna_seq), i) for i in range(0, 3)]
    return res_orig + res_rc

def all_proteins_rf(aa_seq):
    """ Computes all possible proteins in an amino acid sequence. Returns list of possible proteins. """
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins

def all_orfs(dna_seq):
    """ Computes all possible proteins for all open reading frames."""
    rfs = reading_frames(dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots: 
            res.append(p)
    return res

if __name__ == "__main__":
    helloworld()
    import sys
    if len(sys.argv) > 1: 
        in_args= sys.argv[1:] # 첫번째를 제외하고 이후부터 사용
        inFilename = in_args[0]
        outFilename = in_args[1]
        fileCapitalize(inFilename, outFilename)
    else:
        fileCapitalize('myfile2.txt','myfile2_cap.txt')
