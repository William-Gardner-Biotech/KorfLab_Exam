import argparse

# K-mer locations

parser = argparse.ArgumentParser(description='K-mer location program')

parser.add_argument('--file', '-f', required=True, type=str, default=None,
    metavar='<str>', help='Fasta file to be analyzed')

parser.add_argument('--kmer', '-k', required=True, type=int, default=4,
    metavar='<float>', help='K-mer size option')

parser.add_argument('--outfile', '-o', required=False, type=str, default='out_kmers.txt',
    metavar='<str>', help='Name of output file')

parser.add_argument('--negative', '-n', required=False, type=bool, default=False,
    metavar='<bool>', help="Option to count both directions of the strand. [DEFAULT = False]")

arg = parser.parse_args()


with open(arg.file, 'r') as f:
    f.readline()
    sequence = f.read().strip()
    #print(sequence)

k_dict = {}
iter = 0
for nuc in sequence:
    curmer = sequence[iter:iter+arg.kmer]
    if curmer in k_dict:
        k_dict[curmer].append(iter+1)
    else:
        k_dict[curmer] = [(iter+1)]
    if iter == len(sequence)-arg.kmer:
        break
    iter += 1

# Working, on the negative you begin reading it backwards so using conventional
# Right to left logic as you count to the right just make it the negative
if arg.negative:
    right_point = len(sequence)
    # work backwards using the range for loop and subtracting
    for i in range(len(sequence)):
        #print(i)
        left_point = (right_point-arg.kmer)
        negative_mer = (sequence[left_point:right_point])
        # Flips the sub_string
        negative_mer = negative_mer[::-1]
        #print(right_point, negative_mer)
        if negative_mer in k_dict:
            k_dict[negative_mer].append((-1*right_point))
        else:
            k_dict[negative_mer] = [(-1*right_point)]
        if i == len(sequence) - arg.kmer:
            break
        right_point -= 1


output = ''
for key in k_dict.keys():
    output += f'{key} {k_dict[key]}\n'

print(output)
with open(arg.outfile, 'w') as out:
    out.write(output)
