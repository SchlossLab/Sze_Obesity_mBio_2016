$ nano test

#Create new directories for chimera completed files
mkdir chimera.completed
mkdir mothurLogfiles

#Make loop to run the chimera step and move completed samples to respective directory
for chim in $(ls *.fasta)
do
	CFile=${chim//fasta/count_table}
	/share/scratch/schloss/mothur/mothur "#chimera.uchime(fasta=$chim, count=$CFile, dereplicate=t, processors=8)"
	
	mv $chim chimera.completed
	mv $CFile chimera.completed
	mv *logfile mothurLogfiles
	
done

exit 0