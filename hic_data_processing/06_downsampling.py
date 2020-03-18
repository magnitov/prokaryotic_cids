import argparse
import random
import numpy as np

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('input_contacts', type = str, help = 'File with input contacts')
parser.add_argument('number_of_contacts', type = int, help = 'Number of contacts to subsample')
parser.add_argument('output_contacts', type = str, help = 'File with output contacts')
args = parser.parse_args()

input_contacts = args.input_contacts
N = args.number_of_contacts
output_contacts = args.output_contacts

# Read list of contact bins
print('Step 1/5, reading input...')
with open(input_contacts, 'r') as f:
    input_lines = f.readlines()

# Sample so that 1 element represents one contact
print('Step 2/5, splitting input...')
out_lines = []
for l in input_lines:
    num = int(l.split()[-1].replace('\n', ''))
    for j in range(0, num):
        out_lines.append(' '.join(l.split()[:-1]) + ' 1\n')

# Shuffle the final list
print('Step 3/5, shuffling contacts...')
for i in range(5):
    random.shuffle(out_lines)

# Subsample required contacts
print('Step 4/5, subsampling...')
sampled_lines = random.sample(out_lines, N)

# Merge to bins
print('Step 5/5, re-binning contacts...')
sampled_lines_unique = np.unique(sampled_lines, return_counts = True)
write_lines = [str(' '.join(x.split()[:-1])) + ' ' + str(y) + '\n' for (x, y) in zip(sampled_lines_unique[0], sampled_lines_unique[1])]

# Save output
with open(output_contacts, 'w') as g:
    for l in write_lines:
        g.write(l)
