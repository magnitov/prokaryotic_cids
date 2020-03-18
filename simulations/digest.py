import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO
from Bio.Restriction import *
Entrez.email = ""

data = pd.read_csv('prokaryotes.csv')

# Get taxonomic data
data['Taxon'] = [x.split(';')[0] for x in data['Organism Groups'].values]
data['Group'] = [x.split(';')[-1] for x in data['Organism Groups'].values]

# Get chromosomes refseq IDs
chroms = []
delete = []
i = 0
for l in [x.split(';') for x in data['Replicons'].values]:
    new_vals = []
    for val in l:
        if 'chromosome' in val:
            new_vals.append(val.split(':')[1].split('/')[0])
    if len(new_vals) == 1:
        chroms.append(new_vals[0])
    else:
        delete.append(i)
    i += 1
data = data.drop(delete, axis = 0)
data['Chromosome'] = chroms

data = data.drop(['Organism Groups', 'Replicons'], axis = 1)
data.columns = ['Organism', 'Size', 'GC', 'Taxon', 'Group', 'Chromosome']
data = data[['Organism',  'Taxon', 'Group', 'Chromosome', 'Size', 'GC']]
    
def digest(genome_id, restriction_enzyme):
    print('Loading ' + genome_id + '...')
    handle = Entrez.efetch(db = 'nucleotide', id = genome_id, rettype = 'fasta', retmode = 'text')
    record = SeqIO.read(handle, 'fasta')
    seq = record.seq

    print('Digesting ' + genome_id + '...')
    catalysed = []
    if restriction_enzyme == 'HpaII':
        tmp = HpaII.catalyse(seq)
        for frag in tmp:
            catalysed.append(frag)
    if restriction_enzyme == 'HindIII':
        tmp = HindIII.catalyse(seq)
        for frag in tmp:
            catalysed.append(frag)
    if restriction_enzyme == 'MseI':
        tmp = MseI.catalyse(seq)
        for frag in tmp:
            catalysed.append(frag)
    if restriction_enzyme == 'NcoI':
        tmp = NcoI.catalyse(seq)
        for frag in tmp:
            catalysed.append(frag)

    if len(catalysed)!=0:
        catalysed = [len(x) for x in catalysed]
        out = [np.median(catalysed), np.percentile(catalysed, 25), np.percentile(catalysed, 75)]
    else:
        out = [np.nan, np.nan, np.nan]
    return(out)

# Digest
hpaII_data = []
for i in data['Chromosome'].values:
    print(list(data['Chromosome'].values).index(i)+1, '/', len(data), 'HpaII')
    hpaII_data.append(digest(i, 'HpaII'))
pd.DataFrame(hpaII_data).to_csv('hpaII_data.csv', header = 0, index = 0, sep = '\t')

hindIII_data = []
for i in data['Chromosome'].values:
    print(list(data['Chromosome'].values).index(i)+1, '/', len(data), 'HindIII')
    hindIII_data.append(digest(i, 'HindIII'))
pd.DataFrame(hindIII_data).to_csv('hindIII_data.csv', header = 0, index = 0, sep = '\t')
          
ncoI_data = []
for i in data['Chromosome'].values:
    print(list(data['Chromosome'].values).index(i)+1, '/', len(data), 'NcoI')
    ncoI_data.append(digest(i, 'NcoI'))
pd.DataFrame(ncoI_data).to_csv('ncoI_data.csv', header = 0, index = 0, sep = '\t')
          
mseI_data = []
for i in data['Chromosome'].values:
    print(list(data['Chromosome'].values).index(i)+1, '/', len(data), 'MseI')
    mseI_data.append(digest(i, 'MseI'))
pd.DataFrame(mseI_data).to_csv('mseI_data.csv', header = 0, index = 0, sep = '\t')
