#!/bin/sh

#Create new directories for chimera completed files
#mkdir precluster.completed
#mkdir map.completed
#mkdir mothurLogfiles

#Make loop to run the chimera step and move completed samples to respective directory
#for Pclust in $(ls *.fasta)
#do
#	tempfile=${Pclust//fasta/count_table}
#	CFile=${tempfile//filter.unique/filter}
		
#	/share/scratch/schloss/mothur/mothur "#pre.cluster(fasta=$Pclust, count=$CFile, diffs=1)"
	
#	mv $Pclust precluster.completed
#	mv $CFile precluster.completed
#	mv *logfile mothurLogfiles
#	mv *map map.completed
	
#done


#Create new directories for chimera completed files
#mkdir chimera.completed

#Make loop to run the chimera step and move completed samples to respective directory
for chim in $(ls *.fasta)
do
	CFile=${chim//fasta/count_table}
	
	
	/share/scratch/schloss/mothur/mothur "#chimera.uchime(fasta=$chim, count=$CFile, dereplicate=t)"
	
	mv $chim chimera.completed
	mv $CFile chimera.completed
	mv *logfile mothurLogfiles
	
done

exit 0