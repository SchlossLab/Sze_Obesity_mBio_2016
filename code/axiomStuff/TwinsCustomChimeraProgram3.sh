#!/bin/sh

#Create new directories for chimera completed files
mv *.chimeras chimera.completed
mkdir fasta.completed
mkdir accnos.completed

#Make loop to run the chimera step and move completed samples to respective directory
for chim in $(ls *.fasta)
do
	AccFile=${chim//fasta/denovo.uchime.accnos}

	/share/scratch/schloss/mothur/mothur "#remove.seqs(fasta=$chim, accnos=$AccFile)"
	
	mv $chim fasta.completed
	mv $AccFile accnos.completed
	mv *logfile mothurLogfiles
	
done

exit 0
