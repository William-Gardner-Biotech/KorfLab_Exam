import argparse
import math

# Adds CLI options

parser = argparse.ArgumentParser(description='Entropy filter to identify nucleotides with high entropy values, and mask those with low values')
# Revert the required = True && dealute = None after finished && mask deafult
parser.add_argument('--file', '-f', required=True, type=str, default=None,
    metavar='<str>', help='Fasta file to be entropy_masked, can be single fasta sequence or multi-fasta')

parser.add_argument('--entropy', '-e', required=False, type=float, default=1.4,
    metavar='<float>', help='Modify entropy threshold [DEFAULT = 1.4 bits]')

parser.add_argument('--window', '-w', required=False, type=int, default=11,
    metavar='<int>', help='Window size for calculating entropy and iterating over sequence(s)')

parser.add_argument('--mask', '-m', required=False, action='store_true',
    help="Annotate high entropy regions by changing type case, DEFAULT = replaced with 'N'\
        \n USE: python entropy_filter.py -f your.fasta -m [This will make your output mask by changing case, using program without specifying -m will replace low entropy with 'N']")

parser.add_argument('--outfile', '-o', required=False, type=str, default='output.fasta',
    metavar='<str>', help='Name of output file')

parser.add_argument('--dev', '-d', required=False, action='store_true',
    help="Developer mode to visualize more")

arg = parser.parse_args()

# works as a replacement to SeqIO.parse
# this parser will create a list of tuples, (name, sequence)
def fasta_parser(filename):
    sequence_list = []
    name = ''
    file = open(filename, 'r')
    for line in file:
        if line.startswith(">"):
            if not name:
                name = line
                seq = ''
            else:
                # This will make multifasta handle subsequenct sequences
                yield name, seq
                seq = ''
                name = line
        else:
            seq+=line.strip()
    # Catches the append for a single fasta and the final sequence of a multifasta
    yield name, seq


def prob_dist(subseq):
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nucleotide in subseq:
        if nucleotide in counts:
            counts[nucleotide] += 1
        else:
            return False
    return counts

def entropy_mask(seq, window, entropy, mask, dev):
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
            # 'I(E)=log(base2)(1/pE)'
            I_E = math.log(1/(pE), 2)
            if dev:
                print('Window',subseq, 'Information', round(I_E, 1), 'Probability', round(pE, 1), 'Nucleotide', nucleotide)
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

# This is where it all gets compiled
# Takes the tuples list and adds name first then sequence after
    

output = f""
for name, sequence in fasta_parser(arg.file):
    output += f'{name}{entropy_mask(sequence, arg.window, arg.entropy, arg.mask, arg.dev)}\n'

outfile = open(arg.outfile, 'w')
outfile.write(output)

# Notes on current issues and things to resolve 18Apr2023
# 1. Only handles fasta files with nucleotides in upper case (easy fix)
# 2. Output files are missing the description...possible resolve would be to also assign description variable when parsing then call back with index \
# This also would simplify the program removing the need for a tuple list and instead indexing back on original fasta_seqs list
# 3. Finishes every output with a '\n'