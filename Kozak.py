import argparse
import math
import re
from Bio import SeqIO
import json

# Adds CLI options

parser = argparse.ArgumentParser(description='Kozak consensus sequence generator for E.coli using a reference genome and Genbank Annotation File\
                                 \nThe output is a Position Weight Matrix stored as a json file')
# Revert the required = True && dealute = None after finished && mask deafult
parser.add_argument('--file', '-f', required=True, type=str, default=None,
    metavar='<str>', help='GenBank annotation file to be analyzed [.gbff file]')

parser.add_argument('--refgenome', '-g', required=False, type=str, default='data/GCF_000005845.2_ASM584v2_genomic.fna',
    metavar='<str>', help='Location of your complete E.coli genome .fasta file')

parser.add_argument('--outfile', '-o', required=False, type=str, default='output.json',
    metavar='<str>', help='Name of output file')

parser.add_argument('--upstream', '-u', required=False, type=int, default=9,
    metavar='<int>', help='Length of Upstream region from gene start site to analyze\
        \nNOTE: the gene start is added as well so -u 5 will return a PWM of 8')

arg = parser.parse_args()

with open(arg.refgenome, "r") as handle:
    record = SeqIO.read(handle, "fasta")

# Get the nucleotide sequence as a string and remove whitespace and newlines
# Also use the Seq object function to get the reverse complement 
E_coli_reverse = record.seq.reverse_complement()
E_coli_reverse = str(E_coli_reverse).replace(" ", "").replace("\n", "")
E_coli = str(record.seq).replace(" ", "").replace("\n", "")

'''print('len', len(E_coli))
print(E_coli[0:30])
print(E_coli_reverse[4641622:])'''

# Grabs the positions, important to note that GBFF is indexed from 1 not 0 like Python
def collect_pos(gff):
    position_list = []
    features = open(gff, 'r')
    for line in features:
        if 'CDS' in line: 
            # Must find the complement genes because they require special handling
            regmatch = re.search(r"(?:complement\()?(\d+)\.\.", line)
            if regmatch: 
                pos = int(regmatch.group(1))
                # * (-1) is an easy way to mark them as complements
                if 'complement' in line: position_list.append(-1*pos)
                else: position_list.append(pos)
        else:
            continue
    return position_list

# How many bases upstream of start codon should I use, I will do 20 upstream and maybe add a feature to trim it down

def seq_grab(positions, genome, rev_genome, upstream_length):
    kozaks = []
    # genome length is important for deriving start of complement genes
    genome_L = len(genome)
    for pos in positions:
        # handling the positives and negatives
        if pos > 0:
            # skip exception to prevent rare case of gene that starts right next to genome beginning and will cause issues to PWM
            if pos < upstream_length:
                continue
            # pos + 2 is used to grab the start codon (should be ~ATG)
            kozaks.append(genome[pos-1-upstream_length:pos+2])
        else:
            # Opposite of skip for the positives
            if pos - upstream_length < (-1*(genome_L)):
                continue
            else:
                # pos + 1 because GBFF pos is indexed off 1, also pos will be negative
                # We have three variables adding 1 to the index, 1 for GBFF pos, 1 for len(), and 1 for upstream length
                comp_start = genome_L+pos
                kozaks.append(rev_genome[comp_start-upstream_length-1:comp_start+2])
    return kozaks

def PWM(kozaks, upstream_length):
    # init the dictionary of dictionaries, O(n) n^2
    # Outer dict is for the positions and inner is for the nucleotide count probabilities
    # Upstream length assumes counting from 1 not 0 so we must account for that
    outer_dict = {}
    for position in range(0, upstream_length+3):
        up_one = str(position+1)
        inner_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
        replacement_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
        # We must subtract one for the proper indexing
        for kzk in kozaks:
            nuc = kzk[position]
            inner_dict[nuc] += 1
        for key in inner_dict.keys():
            #print(position, key, inner_dict[key])
            replacement_dict[key] = probability_replace(inner_dict[key], inner_dict)
        outer_dict[up_one] = replacement_dict
    return outer_dict

# Information score is 2-I(E)
# Does passing dicts back and forth through this function slow it, I could just sum(divt.values()) in PWM and pass an int
def probability_replace(indv_val, count_dict):
    # Handles divide by zero error and returns what it would logically output of 0
    if indv_val == 0:
        return 0
    total = sum(count_dict.values())
    pE = indv_val/total
    I_E = math.log(1/(pE), 2)
    # subtract the value by 2 which is max information of 4 outcomes
    pwm_score = 2-I_E
    return pwm_score

# Could write into a main() function but why?
x = collect_pos(arg.file)

y = seq_grab(x, E_coli, E_coli_reverse, arg.upstream)

z = PWM(y, arg.upstream)

with open(arg.outfile, 'w') as out:
    json.dump(z, out, indent=4)

#### NOTES #####
# First note that the default Z sequence is 5’(gcc)gccRccAUGG 3’
# GBFF is indexed from 1
# To push the upstream on further up you must: remove one (-1) from comp_start & Minus (1) from positive pos 

### Dead code ####
# This is used to visualize the nested dictionary
'''    # outer key would be the position
    for key in outer_dict.keys():
        print(key)
        # nkey is the inner key AKA the nucleotide
        for nkey in outer_dict[key]:
            print(nkey, outer_dict[key][nkey])'''