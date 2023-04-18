import argparse
import math
from Bio import SeqIO

# Adds CLI options

parser = argparse.ArgumentParser(description='Entropy filter to identify nucleotides with high entropy values')
# Revert the required = True && dealute = None after finished && mask deafult
parser.add_argument('--file', '-f', required=True, type=str, default=None,
    metavar='<str>', help='Fasta file to be analyzed, can be single fasta sequence or multi-fasta')

parser.add_argument('--entropy', '-e', required=False, type=float, default=1.4,
    metavar='<float>', help='Modify entropy threshold [DEFAULT = 1.4 bits]')

parser.add_argument('--window', '-w', required=False, type=int, default=11,
    metavar='<int>', help='Window size for calculating entropy and iterating over sequence(s)')

parser.add_argument('--mask', '-m', required=False, type=bool, default=False,
    metavar='<bool>', help="Annotate high entropy regions by changing type case, DEFAULT = replaced with 'N'")

parser.add_argument('--outfile', '-o', required=False, type=str, default='output.fasta',
    metavar='<str>', help='Name of output file')

arg = parser.parse_args()

fasta_seqs = SeqIO.parse(open(arg.file), 'fasta')
# Builds a list of tuples of all the sequences, helps with multi-fasta, and will maintain correct order
seq_list = []
for fasta in fasta_seqs:
    name, sequence = fasta.id, str(fasta.seq)
    seq_list.append((name, sequence))

# Now we'll go through the sequences and replace the values with the newly updated values

def prob_dist(subseq):
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nucleotide in subseq:
        if nucleotide in counts:
            counts[nucleotide] += 1
        else:
            return False
    return counts

def analyze(seq, window, entropy, mask):
    iter = 0
    new_seq = ''
    for nucleotide in seq:
        subseq = seq[iter:(iter+window)]
        probs = prob_dist(subseq)
        # error handling of non-nucleotide in fasta sequence
        if not probs:
            exit('Improper nucleotide detected while processing')
        else:
            pE = probs[nucleotide]/len(subseq)
            pee = 1/pE
            # 'I(E)=log(base2)(1/pE)'
            I_E = math.log(1/(pE), 2)
            print(subseq, I_E, pee, pE, nucleotide)
            if I_E < entropy:
                if mask:
                    new_seq+=nucleotide.lower()
                else:
                    new_seq+='N'
            else:
                new_seq+=nucleotide
        # Simple way to make final window stop increasing and self resolve 
        if iter >= len(seq)-window:
            pass
        else:
            iter += 1
    return new_seq



analyze(seq_list[0][1], arg.window, arg.entropy, arg.mask)

output = f""
for i in seq_list:
    output += f'{i[0]}\n'
    output += f'{analyze(i[1], arg.window, arg.entropy, arg.mask)}\n'


outfile = open(arg.outfile, 'w')
outfile.write(output)