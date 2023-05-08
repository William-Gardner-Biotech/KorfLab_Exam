import argparse
import math
import re
from Bio import SeqIO

# Adds CLI options

parser = argparse.ArgumentParser(description='Kozak consensus sequence generator for E.coli using a reference genome and Genbank Annotation File')
# Revert the required = True && dealute = None after finished && mask deafult
parser.add_argument('--file', '-f', required=True, type=str, default=None,
    metavar='<str>', help='GenBank annotation file to be analyzed [.gbff file]')

parser.add_argument('--refgenome', '-g', required=False, type=str, default='data/GCF_000005845.2_ASM584v2_genomic.fna',
    metavar='<str>', help='Location of your complete E.coli genome .fasta file')

parser.add_argument('--outfile', '-o', required=False, type=str, default='output.fasta',
    metavar='<str>', help='Name of output file')

parser.add_argument('--length', '-l', required=False, type=int, default=22,
    metavar='<int>', help='Length of ')

parser.add_argument('--dev', '-d', required=False, type=bool, default=False,
    metavar='<bool>', help="Developer mode to visualize more")

arg = parser.parse_args()

with open(arg.refgenome, "r") as handle:
    record = SeqIO.read(handle, "fasta")

# Get the nucleotide sequence as a string and remove whitespace and newlines
E_coli = str(record.seq).replace(" ", "").replace("\n", "")

'''with open(arg.refgenome, 'r') as ref:
    id = ref.readline()
    E_coli = ref.read()
    E_coli.strip('\n')
    E_coli.strip(' ')'''

#print(E_coli)

def collect_pos(gff):
    position_list = []
    features = open(gff, 'r')
    for line in features:
        if 'CDS' in line: 
            regmatch = re.search(r"(?:complement\()?(\d+)\.\.", line)
            if regmatch: 
                pos = int(regmatch.group(1))
                if 'complement' in line: position_list.append(-1*pos)
                else: position_list.append(pos)
        else:
            continue
    return position_list


def seq_grab(positions, genome, length):
    for pos in positions:
        if pos > 0:
            print(genome[pos-1:pos+2])

x = collect_pos(arg.file)
print(x)

seq_grab(x, E_coli, arg.length)


# First note that the default Z sequence is 5’(gcc)gccRccAUGG 3’
# GBFF is indexed from 1