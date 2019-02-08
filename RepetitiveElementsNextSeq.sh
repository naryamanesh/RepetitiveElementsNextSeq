# RepetitiveElementsNextSeq.sh
# The pipeline for finding genomic regions with sequencing problem in NextSeq platform for repetitive elements in whole genome sequencing or targeted sequencing data.
# Created by Nader Aryamanesh and Leila Eshraghi
# Correspondings: Nader.Aryamanesh@uwa.edu.au
# Date uploaded: 08/02/2019

#######################################################################################################################
# Step1: Go to the working directory where the fastq files downloaded.
########################################################################################################################
cd /PATH_TO/rawdata

########################################################################################################################
# Step2: NextSeq platform generates 8 files for each sample so combine the NextSeq files into two files.
########################################################################################################################
ls *L001_R1_001.fastq.gz | sed 's/L001_R1_001.fastq.gz//' > list.txt
cat list.txt | while read line
do
	echo $line
	cat $line$"L00"?"_R1_001.fastq.gz" > $line$"R1.fastq.gz"
	cat $line$"L00"?"_R2_001.fastq.gz" > $line$"R2.fastq.gz"
done

# NextSeq platform has a tendency to insert multiple "G"s to the end of reads with low quality.
#for example:
# @NB551004:340:H2Y75AFXY:1:11101:9558:3809 1:N:0:TCCGGAGA+GCCTCTAT
# CACAGGGTAAACCACCGCCTATCAGGCCCCTGACTGATTCTACCATAGCGGCCAACGATAGACAAGAGCTCGGGCGGGAGGGGGGGGGGCGGCGGGCGGGGGGGGGGGGGCGGCGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGG
# +
# AAAAAEEEEE6EEEEEEEEEEEEEEEEEEAAEEAEE6<EEE//<A//E/EEEE//EE/E6E/EE6EE//////////E//////E/E///////A/////EA///E/E/E////E/<E/E/A//</<<///E/A<////////EE/EA/EE

########################################################################################################################
# Step3: pull out all reads with a minimum of 8 Gs in their sequence.
########################################################################################################################
cat list.txt | while read line
do
	echo $line
	zgrep -B 1 -A 2 "GGGGGGGG" $line$"R1.fastq.gz" | grep -v "^--" > $line$"R1.G8.fastq"
	zgrep -B 1 -A 2 "GGGGGGGG" $line$"R2.fastq.gz" | grep -v "^--" > $line$"R2.G8.fastq"
done

########################################################################################################################
# Step4: quality control of pulled out reads which are "G" rich
########################################################################################################################
mkdir fastqc
module load fastqc/0.11.8 # or any other version
fastqc -o /wrk/aryamane/DONOTREMOVE/Nader_WGS/Nader_WGS-NextSeq/reads/fastqc/ --noextract -f fastq *.G8.fastq

########################################################################################################################
# Step5: convert fastq files into fasta files
########################################################################################################################
cat list.txt | while read line
do
	echo $line
	grep -A 1 "@" $line$"R1.G8.fastq" | sed 's/@/>@/' | grep -v "^--" > $line$"R1.G8.fasta"
	grep -A 1 "@" $line$"R2.G8.fastq" | sed 's/@/>@/' | grep -v "^--" > $line$"R2.G8.fasta"
	cat $line$"R1.G8.fasta" $line$"R2.G8.fasta" > $line$"G8.fasta" # you can combine the R1 and R2 reads together.
done

########################################################################################################################
# Step6: using blastn function to find the reads on the exome/genes. If you want to work genome-wide, you may skip this Step.
########################################################################################################################
# a- create database for exome/genes (in this case fatsa file for Arabidopsis lyrata transcriptome data)
module load biokit
module load emboss
seqret transcriptome_Alyrata.fa transcriptome_Alyrata_db.fasta -osf ncbi
makeblastdb -in transcriptome_Alyrata_db.fasta -out transcriptome_Alyrata_db_ncbi -parse_seqids -dbtype nucl

# b- run blast for files against A. lyrata genes.
cat list.txt | while read line
do
	echo $line
	blastn -task blastn -query $line$"R1.G8.fasta" -db transcriptome_Alyrata_db_ncbi -out $line$"R1.G8.OF6.txt" -outfmt 6 -evalue 5e-2 -perc_identity 80 -max_target_seqs 1 -word_size 40
	blastn -task blastn -query $line$"R2.G8.fasta" -db transcriptome_Alyrata_db_ncbi -out $line$"R2.G8.OF6.txt" -outfmt 6 -evalue 5e-2 -perc_identity 80 -max_target_seqs 1 -word_size 40
	blastn -task blastn -query $line$"G8.fasta" -db transcriptome_Alyrata_db_ncbi -out $line$"G8.OF6.txt" -outfmt 6 -evalue 5e-2 -perc_identity 80 -max_target_seqs 1 -word_size 40
done

# BLASTn tabular output format 6
# Column headers:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# example output:
# @NB551004:340:H2Y75AFXY:1:11101:11267:1223      unk:chr2        97.403  154     1       1       1       151     9168570 9168417 2.31e-68        259
# @NB551004:340:H2Y75AFXY:1:11101:11267:1223      unk:chr2        96.178  157     0       2       1       151     8768892 8769048 9.84e-67        253

########################################################################################################################
# Step7: blast function to find the reads on the chromosomes.
########################################################################################################################
# a- create database for Arabidopsis lyrata chromosomes.
module load biokit
module load emboss
seqret Al_Chrs.fa Al_Chrs_db.fasta -osf ncbi
makeblastdb -in Al_Chrs_db.fasta -out Al_Chrs_db_ncbi -parse_seqids -dbtype nucl

# b- run blast for files against A. lyrata chromosomes
cat list.txt | while read line
do
	echo $line
	blastn -task blastn -query $line$"R1.G8.fasta" -db Al_Chrs_db_ncbi -out $line$"R1.G8.Chrs.OF6.txt" -outfmt 6 -evalue 5e-2 -perc_identity 80 -max_target_seqs 1 -word_size 40
	blastn -task blastn -query $line$"R2.G8.fasta" -db Al_Chrs_db_ncbi -out $line$"R2.G8.Chrs.OF6.txt" -outfmt 6 -evalue 5e-2 -perc_identity 80 -max_target_seqs 1 -word_size 40
	blastn -task blastn -query $line$"G8.fasta" -db Al_Chrs_db_ncbi -out $line$"G8.Chrs.OF6.txt" -outfmt 6 -evalue 5e-2 -perc_identity 80 -max_target_seqs 1 -word_size 40
done

########################################################################################################################
# Step8: convert output files to bed format
########################################################################################################################
cat list.txt | while read line
do
	echo $line
	awk  '{FS="\t"; OFS="\t"; if ($9<$10)  {print $2,$9,$10,$1,$3,$4,$5,$6,$7,$8,$11,$12} else {print $2,$10,$9,$1,$3,$4,$5,$6,$7,$8,$11,$12}}' < $line$"R1.G8.Chrs.OF6.txt" | sed 's/unk:chr//'  > $line$"R1.G8.Chrs.OF6.bed"
	awk  '{FS="\t"; OFS="\t"; if ($9<$10)  {print $2,$9,$10,$1,$3,$4,$5,$6,$7,$8,$11,$12} else {print $2,$10,$9,$1,$3,$4,$5,$6,$7,$8,$11,$12}}' < $line$"R2.G8.Chrs.OF6.txt" | sed 's/unk:chr//'  > $line$"R2.G8.Chrs.OF6.bed"
	awk  '{FS="\t"; OFS="\t"; if ($9<$10)  {print $2,$9,$10,$1,$3,$4,$5,$6,$7,$8,$11,$12} else {print $2,$10,$9,$1,$3,$4,$5,$6,$7,$8,$11,$12}}' < $line$"G8.Chrs.OF6.txt" | sed 's/unk:chr//'  > $line$"G8.Chrs.OF6.bed"
done

########################################################################################################################
# Step9: creating fasta files from bed files, and excluding the sequences with real "G"s in the reference genome
########################################################################################################################
module load gcc/4.7.2
module load bedtools/2.17.0
cat list.txt | while read line
do
	echo $line
	# For R1 reads
	bedtools sort -i $line$"R1.G8.Chrs.OF6.bed" > $line$"R1.G8.Chrs.OF6.sorted.bed"
	bedtools getfasta -fi Al_Chrs.fa -bed $line$"R1.G8.Chrs.OF6.sorted.bed"  -fo $line$"R1.G8.Chrs.OF6.sorted.fa"
	# find the reads with Gs in the reference genome (these reads are good quality and they should be seperated from bad quality reads) 
	grep -B 1 "GGGGGGGG" $line$"R1.G8.Chrs.OF6.sorted.fa"  > $line$"R1.G8.Chrs.OF6.sorted.8G.fa"
	grep -B 1 "CCCCCCCC" $line$"R1.G8.Chrs.OF6.sorted.fa"  > $line$"R1.G8.Chrs.OF6.sorted.8C.fa"
	cat $line$"R1.G8.Chrs.OF6.sorted.8G.fa" $line$"R1.G8.Chrs.OF6.sorted.8C.fa" > $line$"R1.G8.Chrs.OF6.sorted.8GC.fa"
	rm -f $line$"R1.G8.Chrs.OF6.sorted.8G.fa" $line$"R1.G8.Chrs.OF6.sorted.8C.fa"
	grep ">" $line$"R1.G8.Chrs.OF6.sorted.8GC.fa" | sed 's/>//' | sed 's/:/\t/' | sed 's/-/\t/' > $line$"R1.G8.Chrs.OF6.sorted.8GC.bed"
	bedtools merge -i $line$"R1.G8.Chrs.OF6.sorted.8GC.bed" > $line$"R1.G8.Chrs.OF6.sorted.8GC.merged.bed"
	# subtract the good quality read coordinates from total reads with Gs.
	bedtools subtract -a $line$"R1.G8.Chrs.OF6.sorted.bed" -b $line$"R1.G8.Chrs.OF6.sorted.8GC.merged.bed" -A > $line$"R1.G8.Chrs.OF6.sorted.sub8GC.bed"
	# merge regions with 200bp gap as the libraries sizes were around 400bp.
	bedtools merge -d 200 -n -i  $line$"R1.G8.Chrs.OF6.sorted.sub8GC.bed" > $line$"R1.G8.Chrs.OF6.sorted.sub8GC.merged.bed"
	bedtools getfasta -fi Al_Chrs.fa -bed $line$"R1.G8.Chrs.OF6.sorted.sub8GC.merged.bed"  -fo $line$"R1.G8.Chrs.OF6.sorted.sub8GC.merged.fa"

	# For R2 reads
	bedtools sort -i $line$"R2.G8.Chrs.OF6.bed" > $line$"R2.G8.Chrs.OF6.sorted.bed"
	bedtools getfasta -fi Al_Chrs.fa -bed $line$"R2.G8.Chrs.OF6.sorted.bed"  -fo $line$"R2.G8.Chrs.OF6.sorted.fa"
	
	grep -B 1 "GGGGGGGG" $line$"R2.G8.Chrs.OF6.sorted.fa"  > $line$"R2.G8.Chrs.OF6.sorted.8G.fa"
	grep -B 1 "CCCCCCCC" $line$"R2.G8.Chrs.OF6.sorted.fa"  > $line$"R2.G8.Chrs.OF6.sorted.8C.fa"
	cat $line$"R2.G8.Chrs.OF6.sorted.8G.fa" $line$"R2.G8.Chrs.OF6.sorted.8C.fa" > $line$"R2.G8.Chrs.OF6.sorted.8GC.fa"
	rm -f $line$"R2.G8.Chrs.OF6.sorted.8G.fa" $line$"R2.G8.Chrs.OF6.sorted.8C.fa"
	grep ">" $line$"R2.G8.Chrs.OF6.sorted.8GC.fa" | sed 's/>//' | sed 's/:/\t/' | sed 's/-/\t/' > $line$"R2.G8.Chrs.OF6.sorted.8GC.bed"
	bedtools merge -i $line$"R2.G8.Chrs.OF6.sorted.8GC.bed" > $line$"R2.G8.Chrs.OF6.sorted.8GC.merged.bed"
	
	bedtools subtract -a $line$"R2.G8.Chrs.OF6.sorted.bed" -b $line$"R2.G8.Chrs.OF6.sorted.8GC.merged.bed" -A > $line$"R2.G8.Chrs.OF6.sorted.sub8GC.bed"
	bedtools merge -d 200 -n -i  $line$"R2.G8.Chrs.OF6.sorted.sub8GC.bed" > $line$"R2.G8.Chrs.OF6.sorted.sub8GC.merged.bed"
	bedtools getfasta -fi Al_Chrs.fa -bed $line$"R2.G8.Chrs.OF6.sorted.sub8GC.merged.bed"  -fo $line$"R2.G8.Chrs.OF6.sorted.sub8GC.merged.fa"

	# For combined R1 and R2 reads
	bedtools sort -i $line$"G8.Chrs.OF6.bed" > $line$"G8.Chrs.OF6.sorted.bed"
	bedtools getfasta -fi Al_Chrs.fa -bed $line$"G8.Chrs.OF6.sorted.bed"  -fo $line$"G8.Chrs.OF6.sorted.fa"
	
	grep -B 1 "GGGGGGGG" $line$"G8.Chrs.OF6.sorted.fa"  > $line$"G8.Chrs.OF6.sorted.8G.fa"
	grep -B 1 "CCCCCCCC" $line$"G8.Chrs.OF6.sorted.fa"  > $line$"G8.Chrs.OF6.sorted.8C.fa"
	cat $line$"G8.Chrs.OF6.sorted.8G.fa" $line$"G8.Chrs.OF6.sorted.8C.fa" > $line$"G8.Chrs.OF6.sorted.8GC.fa"
	rm -f $line$"G8.Chrs.OF6.sorted.8G.fa" $line$"G8.Chrs.OF6.sorted.8C.fa"
	grep ">" $line$"G8.Chrs.OF6.sorted.8GC.fa" | sed 's/>//' | sed 's/:/\t/' | sed 's/-/\t/' > $line$"G8.Chrs.OF6.sorted.8GC.bed"
	bedtools merge -i $line$"G8.Chrs.OF6.sorted.8GC.bed" > $line$"G8.Chrs.OF6.sorted.8GC.merged.bed"
	
	bedtools subtract -a $line$"G8.Chrs.OF6.sorted.bed" -b $line$"G8.Chrs.OF6.sorted.8GC.merged.bed" -A > $line$"G8.Chrs.OF6.sorted.sub8GC.bed"
	bedtools merge -d 200 -n -i  $line$"G8.Chrs.OF6.sorted.sub8GC.bed" > $line$"G8.Chrs.OF6.sorted.sub8GC.merged.bed"
	bedtools getfasta -fi Al_Chrs.fa -bed $line$"G8.Chrs.OF6.sorted.sub8GC.merged.bed"  -fo $line$"G8.Chrs.OF6.sorted.sub8GC.merged.fa"
done

########################################################################################################################
# Step10: Combine all samples into a single file and get fasta file
########################################################################################################################
cat *_G8.Chrs.OF6.sorted.sub8GC.bed > NextSeq.all.G8.Chrs.OF6.sub8GC.bed
bedtools sort -i NextSeq.all.G8.Chrs.OF6.sub8GC.bed > NextSeq.all.G8.Chrs.OF6.sub8GC.sorted.bed
bedtools merge -d 200 -n -i NextSeq.all.G8.Chrs.OF6.sub8GC.sorted.bed > NextSeq.all.G8.Chrs.OF6.sub8GC.sorted.merged.bed
bedtools getfasta -fi /wrk/aryamane/DONOTREMOVE/iGenomes/Arabidopsis_lyrata/Al_Chrs.fa -bed NextSeq.all.G8.Chrs.OF6.sub8GC.sorted.merged.bed  -fo NextSeq.all.G8.Chrs.OF6.sub8GC.sorted.merged.fa

########################################################################################################################
# Step11: count number of location with repetative elements in each chromosome
########################################################################################################################
cat list.txt | while read line
do
	echo $line >> repetative_count.sub8GC.txt
	for i in {1..8}
	do	
	echo "Number of locations with repetative elements for file R1 on chromosome" ${i}  >> repetative_count.sub8GC.txt
	grep -cw ${i} $line$"R1.G8.Chrs.OF6.sorted.sub8GC.merged.bed" >> repetative_count.sub8GC.txt
	echo "Number of locations with repetative elements for file R2 on chromosome " ${i}  >> repetative_count.sub8GC.txt
	grep -cw ${i} $line$"R2.G8.Chrs.OF6.sorted.sub8GC.merged.bed" >> repetative_count.sub8GC.txt
	echo "Number of locations with repetative elements for combined file on chromosome" ${i}  >> repetative_count.sub8GC.txt
	grep -cw ${i} $line$"G8.Chrs.OF6.sorted.sub8GC.merged.bed" >> repetative_count.sub8GC.txt
	done
done

########################################################################################################################
# Step12: plot repeat density in R
########################################################################################################################
module load r-env
R

library(ggplot2)
# for individual samples:
samples=c("Lib-1_S1_","Lib-2_S2_", "Lib-3_S3_", "Lib-4_S4_", "Lib-5_S5_", "Lib-6_S6_")
for(j in 1:6){
repeats<-read.table(paste(samples[j], "G8.Chrs.OF6.sorted.sub8GC.bed", sep=''),sep="\t",header=F)
head(repeats)
colnames(repeats)<-c("chr","start","end")
head(repeats)

repeatsDensity<-ggplot(repeats) + 
geom_histogram(aes(x=start),binwidth=1e3) + 
facet_wrap(~ chr,ncol=1, scales = "free") +
ggtitle("Density of repeats with sequencing problem across Arabidopsis lyrata chromosomes") +
scale_y_log10() +
xlab("Position in the genome") + 
ylab("repeats density") + 
theme_bw() # I prefer the black and white theme
pdf(file = paste("Repeats density for ",samples[j],"merged_sub8GC_1k.pdf",sep=''),width=10,height=15)
print(repeatsDensity)
dev.off()
}

# for comined sample
repeats<-read.table("NextSeq.all.G8.Chrs.OF6.sub8GC.sorted.bed", sep="\t",header=F)
head(repeats)
colnames(repeats)<-c("chr","start","end")
head(repeats)

repeatsDensity<-ggplot(repeats) + 
geom_histogram(aes(x=start),binwidth=1e3) + 
facet_wrap(~ chr,ncol=1, scales = "free") +
ggtitle("Density of repeats with sequencing problem across Arabidopsis lyrata chromosomes") +
scale_y_log10() +
xlab("Position in the genome") + 
ylab("repeats density") + 
theme_bw() # I prefer the black and white theme
pdf(file = "Repeats density for NextSeq_all_merged_sub8GC_1ksubs.pdf",width=10,height=15)
print(repeatsDensity)
dev.off()
q()
n

########################################################################################################################
# Step13: visualize the fasta file for a region of interest with the repeat finder "REPFIND":
########################################################################################################################

# Example input
# find a coordinate of ineterest in the bed file (e.g. NextSeq.all.G8.Chrs.OF6.sub8GC.sorted.merged.bed)
# 1	21784930	21785057	3

# find the sequence of the coordinate in the fasta file (e.g. NextSeq.all.G8.Chrs.OF6.sub8GC.sorted.merged.fa)
# >1:21784930-21785057
# TACCGAGAGCTCGTAAGGTTATACAGCCCAGAGTGTGAGTGGTATGGTTATATACTCAAAGCTCGTAAGGTTATATAGCTCAGAGTGTGAGTGGTACGGATATGATACCCAGAGCTCGTAAGATTAT

# go to https://zlab.bu.edu/repfind/form.html
# copy and paste the sequence of ineterst in the box (e.g. Minimum repeat length of 6)
# EXAMPLE output:
#Word: AGCTCGTAAG
# Word locations: 8 60 113 
# Most significant cluster: 8 to 122
# P-value: 1.64589e-08
# Word: AGCTCGTAAGGTTATA
# Word locations: 8 60 
# Most significant cluster: 8 to 75
# P-value: 4.39185e-08
# Word: CAGAGTGTGAGTGGTA
# Word locations: 29 81 
# Most significant cluster: 29 to 96
# P-value: 1.94576e-07
# Word: GGTTATA
# Word locations: 17 46 69 
# Most significant cluster: 17 to 75
# P-value: 1.97892e-05
# Word: AGAGCTCGTAAG
# Word locations: 6 111 
# Most significant cluster: 6 to 122
# P-value: 2.24601e-05
# **********************************
