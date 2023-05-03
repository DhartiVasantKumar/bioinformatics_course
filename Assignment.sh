#!/bin/bash #

#The basic tools and packages are installed using conda
cd ~/

wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh

chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh

bash ./Anaconda3-2022.10-Linux-x86_64.sh

source ~/.bashrc

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib
#The data organisation is carried out as follows by creating directories and files for storing data
pwd
mkdir ngs_course
mkdir ngs_course/dnaseq
cd ngs_course/dnaseq
mkdir data meta results logs
ls -lF
cd ~/ngs_course/dnaseq/data
mkdir untrimmed_fastq
mkdir trimmed_fastq

#The input files are downloaded using wget according to the links given

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
mv *fastq.qz ~/ngs_course/dnaseq/data/untrimmed_fastq

#The files are unzipped to carry out quality cotrol and alignment

cd untrimmed_fastq 
zcat NGS0001.R1.fastq.qz > NGS0001.R1.fastq
zcat NGS0001.R2.fastq.qz > NGS0001.R2.fastq
head NGS0001.R1.fastq
head NGS0001.R2.fastq

#Fastqc is carried out on the fastq data

fastqc -t 4 *.fastq.qz
mkdir ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/
ls
cd
cd ~/ngs_course/dnaseq/results/
cd ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/
ls
for zip in *.zip
> do
> unzip $zip
> done
cat */summary.txt > ~/ngs_course/dnaseq/logs/fastqc_summaries.txt
ls -lh NGS0001.R1.fastq.qz_fastqc
ls -lh NGS0001.R2.fastq.qz_fastqc

#Trimmomatic operation is carried out to remove lower quality reads by specifying read group identifiers

cd ~/ngs_course/dnaseq/logs/
cd ~/ngs_course/dnaseq/data/untrimmed_fastq
$ trimmomatic PE -threads 4 -phred33 /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq -baseout /home/ubuntu/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50

# Alignment

## Read alignment using bwa mem while inputting the correct read group identifiers

mkdir -p ~/ngs_course/dnaseq/data/reference
mv ~/ngs_course/dnaseq/data/hg19.fa.gz ~/ngs_course/dnaseq/data/reference/
bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:11V6WR1' -I 250,50  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/trimmed_
fastq/NGS0001_trimmed_R_1P ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/ng

## Convert the file to bam file and generate index using samtools

samtools view -h -b NGS0001.sam > NGS0001.bam
samtools sort NGS0001.bam > NGS0001_sorted.bam
samtools index NGS0001_sorted.bam > NGS0001_sorted.bam.bai

# Post Alignment QC and filtering

## MarkDuplicates

picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_sorted_marked.bam

## Bam file is filtered based on mapping quality and bitwise flags using samtools

samtools index NGS0001_sorted_marked.bam NGS0001_index.bai
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam

## Alignment statistic analysis (flagstats, idxstats)

samtools idxstats NGS0001_sorted.bam
samtools flagstat NGS0001_sorted.bam
samtools stats -i NGS0001_sorted.bam

# Variant calling was carried out using Freebayes

zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa
samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa
freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/results/NGS0001.vcf
bgzip ~/ngs_course/dnaseq/results/NGS0001.vcf
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001.vcf.gz

# Vcf filtering by specifying Freebayes Information fields (low quality calls, low read depth, alleles seen on one strand, etc are removed during filtering)

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
	 > ~/ngs_course/dnaseq/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/resul ts/NGS0001_filtered.vcf
bedtools intersect -header -wa -a ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf -b ~/ngs_course/dnaseq/data/annotation.bed > ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.vcf
bgzip ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.vcf
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.vcf.gz

# Annotating the variants using wANNOVAR

tar -zxvf annovar.latest.tar.gz
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

# Viewing the annotated and filtered variants using Excel via Filezilla

./table_annovar.pl ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.avinput humandb/ -buildver hg19  -out ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

