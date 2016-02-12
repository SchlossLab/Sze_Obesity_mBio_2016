#!/bin/sh

/share/scratch/schloss/mothur/mothur "#merge.files(input=group1.fasta-group2.fasta-group3.fasta-group4.fasta-group10.fasta-group11.fasta-group12.fasta-group13.fasta-group14.fasta-group15.fasta-group16.fasta-group17.fasta-group18.fasta-group19.fasta-group20.fasta-group21.fasta-group22.fasta-group23.fasta-group24.fasta-group25.fasta-group26.fasta-group27.fasta-group28.fasta-group29.fasta-group30.fasta-group31.fasta-group32.fasta-group33.fasta-group34.fasta-group35.fasta-group36.fasta-group37.fasta-group38.fasta-group39.fasta-group40.fasta-group41.fasta-group42.fasta-group43.fasta-group44.fasta-group45.fasta-group46.fasta-group47.fasta-group48.fasta-group49.fasta, output=combined.fasta)"
/share/scratch/schloss/mothur/mothur "#merge.files(input=group1.groups-group2.groups-group3.groups-group4.groups-group10.groups-group11.groups-group12.groups-group13.groups-group14.groups-group15.groups-group16.groups-group17.groups-group18.groups-group19.groups-group20.groups-group21.groups-group22.groups-group23.groups-group24.groups-group25.groups-group26.groups-group27.groups-group28.groups-group29.groups-group30.groups-group31.groups-group32.groups-group33.groups-group34.groups-group35.groups-group36.groups-group37.groups-group38.groups-group39.groups-group40.groups-group41.groups-group42.groups-group43.groups-group44.groups-group45.groups-group46.groups-group47.groups-group48.groups-group49.groups, output=combined.groups)"
/share/scratch/schloss/mothur/mothur "#unique.seqs(fasta=combined.fasta)"
/share/scratch/schloss/mothur/mothur "#count.seqs(name=combined.names, group=combined.groups)"
/share/scratch/schloss/mothur/mothur "#pre.cluster(fasta=combined.unique.fasta, count=combined.count_table, diffs=1)"

mv *logfile mothurLogfiles 


exit 0
