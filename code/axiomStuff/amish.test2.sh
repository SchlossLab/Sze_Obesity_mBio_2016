#!/bin/bash

# Make a directory to store the reference files
#mkdir reference.files
#mv silva.bacteria.fasta reference.files
#mv trainset9_032012.pds.tax reference.files
#mv trainset9_032012.pds.fasta reference.files

# Read in the sra and make a fastq instead of a sff.
#mkdir completed.redo.sra
#mkdir logfile.completed

#for sample in $(ls *sra)
#do
#	fastq-dump $sample
#	mv $sample completed.redo.sra
#	mv *logfile logfile.completed 
	
#done

# Get the quality file and fasta information from the fastq
mkdir completed.fastq

for sample2 in $(ls *fastq)
do
	/share/scratch/schloss/mothur/mothur "#fastq.info(fastq=$sample2)"	
	mv $sample2 completed.fastq
	mv *logfile logfile.completed
	
done

# Trim the fastas according to quality 
mkdir first.trim.completed
for Trimming in $(ls *fasta)
do
	Quality=${Trimming//fasta/qual}
	/share/scratch/schloss/mothur/mothur "#trim.seqs(fasta=$Trimming, qfile=$Quality, maxambig=0, maxhomop=8, qwindowaverage=35, qwindowsize=50)"
	mv $Trimming first.trim.completed
	mv $Quality first.trim.completed
	mv *logfile logfile.completed
done

mkdir scrap.fasta.completed
mkdir scrap.qual.completed
mkdir trim.qual.completed
mv *scrap.qual scrap.qual.completed
mv * scrap.fasta.completed
mv *trim.qual trim.qual.completed


