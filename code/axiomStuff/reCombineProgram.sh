#!/bin/sh

/share/scratch/schloss/mothur/mothur "#merge.files(input=group1.fasta-group2.fasta-group3.fasta-group4.fasta-group5.fasta-group6.fasta-group7.fasta-group8.fasta-group9.fasta-group10.fasta-group11.fasta-group12.fasta-group13.fasta-group14.fasta-group15.fasta-group16.fasta-group17.fasta-group18.fasta-group19.fasta-group20.fasta-group21.fasta-group22.fasta-group23.fasta-group24.fasta-group25.fasta, output=combined.fasta)"
/share/scratch/schloss/mothur/mothur "#merge.files(input=group1.groups-group2.groups-group3.groups-group4.groups-group5.groups-group6.groups-group7.groups-group8.groups-group9.groups-group10.groups-group11.groups-group12.groups-group13.groups-group14.groups-group15.groups-group16.groups-group17.groups-group18.groups-group19.groups-group20.groups-group21.groups-group22.groups-group23.groups-group24.groups-group25.groups, output=combined.groups)"
/share/scratch/schloss/mothur/mothur "#unique.seqs(fasta=combined.fasta)"
/share/scratch/schloss/mothur/mothur "#count.seqs(name=combined.names, group=combined.groups)"

mv *logfile mothurLogfiles 


exit 0
