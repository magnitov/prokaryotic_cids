GENOME=$1
RESOLUTION=$2

# Create folders to store data and genome
mkdir -p bac_STARgenome

# Produce alignment of RNAseq reads
echo 'Start genome preparations...'
/home/magnitov/Software/STAR-2.6.1c/source/./STAR --runMode genomeGenerate --genomeSAindexNbases 8 --runThreadN 16 --genomeDir ./bac_STARgenome --genomeFastaFiles /home/magnitov/data/genomes/${GENOME}.fa
echo 'Start mapping RNA-seq rep1...'
/home/magnitov/Software/STAR-2.6.1c/source/./STAR --genomeDir ./bac_STARgenome --alignIntronMax 1 --runThreadN 24 --outFileNamePrefix  RNAseq_ --outSAMheaderHD @HD --readFilesCommand zcat --readFilesIn RNAseq.fastq.gz

# Manipulate and sort SAM files
echo 'Start manipulations with SAM files...'
samtools view -@ 16 -q 30 -o RNAseq.bam RNAseq_Aligned.out.sam
echo 'Sorting...'
samtools sort -@ 16 -o RNAseq_sorted.bam -T sorting_ -@ 16 RNAseq.bam
echo 'Generating BED intervals...'
/home/magnitov/Software/bedtools2/bin/./bedtools makewindows -g ${GENOME}_chromsizes.txt -w ${RESOLUTION} > bac_windows.bed
/home/magnitov/Software/bedtools2/bin/./bedtools intersect -a bac_windows.bed -b RNAseq_sorted.bam -c -F 0.5 -sorted > RNAseq.bed
echo 'Done!'
