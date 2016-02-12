#!/bin/sh

Create new directories for chimera completed files
mkdir deunique.fasta.completed
mkdir deunique.countTable.completed

Make loop to run the chimera step and move completed samples to respective directory
for chim in $(ls *pick.fasta)
do
	CTable=${chim//pick.fasta/denovo.uchime.pick.count_table}

	/share/scratch/schloss/mothur/mothur "#deunique.seqs(fasta=$chim, count=$CTable)"
	
	mv $chim deunique.fasta.completed
	mv $CTable deunique.countTable.completed
	mv *logfile mothurLogfiles
done
	
mkdir redundant.fasta.completed
mkdir redundant.groups.completed
	
for RFasta in $(ls *pick.redundant.fasta)
do
	RGroups=${RFasta//pick.redundant.fasta/denovo.uchime.pick.redundant.groups}
	
	/share/scratch/schloss/mothur/mothur "#rename.seqs(fasta=$RFasta, group=$RGroups)"

	mv $RFasta redundant.fasta.completed
	mv $RGroups redundant.groups.completed
	mv *renamed_map redundant.groups.completed
	mv *logfile mothurLogfiles
done

export var=1
mkdir redundant.renamed.groups.completed
mkdir redundant.renamed.fasta.completed
echo ${var}
for renamedF in $(ls *pick.redundant.renamed.fasta)
do
	renamedG=${renamedF//pick.redundant.renamed.fasta/denovo.uchime.pick.redundant.renamed.groups}
	echo ${renamedF}
	echo ${renamedG}
	echo ${var}
	
	cp $renamedF ${var}.fasta
	cp $renamedG ${var}.groups
	
	mv $renamedF redundant.renamed.fasta.completed
	mv $renamedG redundant.renamed.groups.completed
	var=$((var+1))
done	
	
	
/share/scratch/schloss/mothur/mothur "#merge.files(input=1.fasta-2.fasta-3.fasta-4.fasta-5.fasta-6.fasta-7.fasta-8.fasta-9.fasta-10.fasta-11.fasta-12.fasta, output=group42.fasta)"
/share/scratch/schloss/mothur/mothur "#merge.files(input=1.groups-2.groups-3.groups-4.groups-5.groups-6.groups-7.groups-8.groups-9.groups-10.groups-11.groups-12.groups, output=group42.groups)"

mv *logfile mothurLogfiles 


exit 0
