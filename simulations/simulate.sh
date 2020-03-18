ENZYME=$1
READ_NUM=30000000

# Bifidobacterium adolescentis
python /home/magnitov/Software/sim3C-0.1.1/sim3C.py -r 1234 -m hic --num-pairs ${READ_NUM} -e ${ENZYME} --efficiency 0.75 --dist uniform -l 50 /home/magnitov/data/genomes/bif_adol1.fa Bifidobacterium_adolescentis_${ENZYME}.fastq
grep -A 3 '1:Y:18:1' Bifidobacterium_adolescentis_${ENZYME}.fastq | grep -v '^--$' > Bifidobacterium_adolescentis_${ENZYME}_R1.fastq
grep -A 3 '2:Y:18:1' Bifidobacterium_adolescentis_${ENZYME}.fastq | grep -v '^--$' > Bifidobacterium_adolescentis_${ENZYME}_R2.fastq
rm profile.tsv Bifidobacterium_adolescentis_${ENZYME}.fastq

# Clostridium difficile
python /home/magnitov/Software/sim3C-0.1.1/sim3C.py -r 5678 -m hic --num-pairs ${READ_NUM} -e ${ENZYME} --efficiency 0.75 --dist uniform -l 50 /home/magnitov/data/genomes/clost_diff1.fa Clostridium_difficile_${ENZYME}.fastq
grep -A 3 '1:Y:18:1' Clostridium_difficile_${ENZYME}.fastq | grep -v '^--$' > Clostridium_difficile_${ENZYME}_R1.fastq
grep -A 3 '2:Y:18:1' Clostridium_difficile_${ENZYME}.fastq | grep -v '^--$' > Clostridium_difficile_${ENZYME}_R2.fastq
rm profile.tsv Clostridium_difficile_${ENZYME}.fastq

# Bacteroides fragilis
python /home/magnitov/Software/sim3C-0.1.1/sim3C.py -r 1589 -m hic --num-pairs ${READ_NUM} -e ${ENZYME} --efficiency 0.75 --dist uniform -l 50 /home/magnitov/data/genomes/bact_frag1.fa Bacteroides_fragilis_${ENZYME}.fastq
grep -A 3 '1:Y:18:1' Bacteroides_fragilis_${ENZYME}.fastq | grep -v '^--$' > Bacteroides_fragilis_${ENZYME}_R1.fastq
grep -A 3 '2:Y:18:1' Bacteroides_fragilis_${ENZYME}.fastq | grep -v '^--$' > Bacteroides_fragilis_${ENZYME}_R2.fastq
rm profile.tsv Bacteroides_fragilis_${ENZYME}.fastq

# Pseudomonas aeruginosa
python /home/magnitov/Software/sim3C-0.1.1/sim3C.py -r 9842 -m hic --num-pairs ${READ_NUM} -e ${ENZYME} --efficiency 0.75 --dist uniform -l 50 /home/magnitov/data/genomes/pseu_aer1.fa Pseudomonas_aeruginosa_${ENZYME}.fastq
grep -A 3 '1:Y:18:1' Pseudomonas_aeruginosa_${ENZYME}.fastq | grep -v '^--$' > Pseudomonas_aeruginosa_${ENZYME}_R1.fastq
grep -A 3 '2:Y:18:1' Pseudomonas_aeruginosa_${ENZYME}.fastq | grep -v '^--$' > Pseudomonas_aeruginosa_${ENZYME}_R2.fastq
rm profile.tsv Pseudomonas_aeruginosa_${ENZYME}.fastq

# Bacillus subtilis
python /home/magnitov/Software/sim3C-0.1.1/sim3C.py -r 7365 -m hic --num-pairs ${READ_NUM} -e ${ENZYME} --efficiency 0.75 --dist uniform -l 50 /home/magnitov/data/genomes/bs1.fa Bacillus_subtilis_${ENZYME}.fastq
grep -A 3 '1:Y:18:1' Bacillus_subtilis_${ENZYME}.fastq | grep -v '^--$' > Bacillus_subtilis_${ENZYME}_R1.fastq
grep -A 3 '2:Y:18:1' Bacillus_subtilis_${ENZYME}.fastq | grep -v '^--$' > Bacillus_subtilis_${ENZYME}_R2.fastq
rm profile.tsv Bacillus_subtilis_${ENZYME}.fastq

# Caulobacter crescentus
python /home/magnitov/Software/sim3C-0.1.1/sim3C.py -r 1083 -m hic --num-pairs ${READ_NUM} -e ${ENZYME} --efficiency 0.75 --dist uniform -l 50 /home/magnitov/data/genomes/cc1.fa Caulobacter_crescentus_${ENZYME}.fastq
grep -A 3 '1:Y:18:1' Caulobacter_crescentus_${ENZYME}.fastq | grep -v '^--$' > Caulobacter_crescentus_${ENZYME}_R1.fastq
grep -A 3 '2:Y:18:1' Caulobacter_crescentus_${ENZYME}.fastq | grep -v '^--$' > Caulobacter_crescentus_${ENZYME}_R2.fastq
rm profile.tsv Caulobacter_crescentus_${ENZYME}.fastq

gzip *.fastq
