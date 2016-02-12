#!/bin/sh

Create new directories for chimera completed files
mkdir uniqueSeqs.completed

Make loop to run the chimera step and move completed samples to respective directory
for unique in $(ls *fasta)
do
	#CTable=${chim//pick.fasta/denovo.uchime.pick.count_table}

	/share/scratch/schloss/mothur/mothur "#unique.seqs(fasta=$unique)"
	
	mv $unique uniqueSeqs.completed
	mv *logfile mothurLogfiles
done
		
/share/scratch/schloss/mothur/mothur "#merge.files(input=group1.unique.fasta-group2.unique.fasta-group3.unique.fasta-group4.unique.fasta-group10.unique.fasta-group11.unique.fasta-group12.unique.fasta-group13.unique.fasta-group14.unique.fasta-group15.unique.fasta-group16.unique.fasta-group17.unique.fasta-group18.unique.fasta-group19.unique.fasta-group20.unique.fasta-group21.unique.fasta-group22.unique.fasta-group23.unique.fasta-group24.unique.fasta-group25.unique.fasta-group26.unique.fasta-group27.unique.fasta-group28.unique.fasta-group29.unique.fasta-group30.unique.fasta-group31.unique.fasta-group32.unique.fasta-group33.unique.fasta-group34.unique.fasta-group35.unique.fasta-group36.unique.fasta-group37.unique.fasta-group38.unique.fasta-group39.unique.fasta-group40.unique.fasta-group41.unique.fasta-group42.unique.fasta-group43.unique.fasta-group44.unique.fasta-group45.unique.fasta-group46.unique.fasta-group47.unique.fasta-group48.unique.fasta-group49.unique.fasta, output=combined.unique.fasta)"
/share/scratch/schloss/mothur/mothur "#merge.files(input=group1.groups-group2.groups-group3.groups-group4.groups-group10.groups-group11.groups-group12.groups-group13.groups-group14.groups-group15.groups-group16.groups-group17.groups-group18.groups-group19.groups-group20.groups-group21.groups-group22.groups-group23.groups-group24.groups-group25.groups-group26.groups-group27.groups-group28.groups-group29.groups-group30.groups-group31.groups-group32.groups-group33.groups-group34.groups-group35.groups-group36.groups-group37.groups-group38.groups-group39.groups-group40.groups-group41.groups-group42.groups-group43.groups-group44.groups-group45.groups-group46.groups-group47.groups-group48.groups-group49.groups, output=combined.groups)"
/share/scratch/schloss/mothur/mothur "#merge.files(input=group1.names-group2.names-group3.names-group4.names-group10.names-group11.names-group12.names-group13.names-group14.names-group15.names-group16.names-group17.names-group18.names-group19.names-group20.names-group21.names-group22.names-group23.names-group24.names-group25.names-group26.names-group27.names-group28.names-group29.names-group30.names-group31.names-group32.names-group33.names-group34.names-group35.names-group36.names-group37.names-group38.names-group39.names-group40.names-group41.names-group42.names-group43.names-group44.names-group45.names-group46.names-group47.names-group48.names-group49.names, output=combined.names)"

mv *logfile mothurLogfiles 

exit 0
