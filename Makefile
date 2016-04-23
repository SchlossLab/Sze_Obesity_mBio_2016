REFS = data/references
FIGS = submission
TABLES = results/tables
PROC = data/process

# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'



################################################################################
#
# Part 1: Get the references
#
# We will need several reference files to complete the analyses including the
# SILVA reference alignment and RDP reference taxonomy.
#
################################################################################

# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v123 and described at http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/
# We will use the SEED v. 123, which contain 12,083 bacterial sequences. This
# also contains the reference taxonomy. We will limit the databases to only
# include bacterial sequences.

$(REFS)/silva.seed.align :
	wget -N http://mothur.org/w/images/1/15/Silva.seed_v123.tgz
	tar xvzf Silva.seed_v123.tgz silva.seed_v123.align silva.seed_v123.tax
	mothur "#get.lineage(fasta=silva.seed_v123.align, taxonomy=silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v123.pick.align, processors=8)"
	mv silva.seed_v123.pick.align $(REFS)/silva.seed.align
	rm Silva.seed_v123.tgz silva.seed_v123.*


# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

$(REFS)/trainset14_032015.% :
	wget -N http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds/trainset14_032015.pds.*
	mv trainset14_032015.pds/* $(REFS)/
	rmdir trainset14_032015.pds
	rm Trainset14_032015.pds.tgz


################################################################################
#
#	Part 2: Get shared file, taxonomy data, and metadata from individual
#	studies
#
################################################################################


# Get Baxter study data from project repository...
data/baxter/baxter.braycurtis.0.03.lt.ave.dist\
	data/baxter/baxter.groups.ave-std.summary\
	data/baxter/baxter.0.03.subsample.shared\
	data/baxter/baxter.metadata\
	data/baxter/baxter.taxonomy : code/baxter.batch code/baxter.R
	bash code/baxter.batch


# Get HMP study data from project repository...
data/hmp/hmp.braycurtis.0.03.lt.ave.dist\
	data/hmp/hmp.groups.ave-std.summary\
	data/hmp/hmp.0.03.subsample.shared\
	data/hmp/hmp.metadata\
	data/hmp/hmp.taxonomy : code/hmp.batch code/hmp.R
	bash code/hmp.batch


# Get Ross study data and process through mothur
data/ross/ross.braycurtis.0.03.lt.ave.dist\
	data/ross/ross.groups.ave-std.summary\
	data/ross/ross.0.03.subsample.shared\
	data/ross/ross.metadata\
	data/ross/ross.taxonomy : code/ross.batch code/ross.R
	bash code/ross.batch


# Get Escobar study data and process through mothur
data/escobar/escobar.braycurtis.0.03.lt.ave.dist\
	data/escobar/escobar.groups.ave-std.summary\
	data/escobar/escobar.0.03.subsample.shared\
	data/escobar/escobar.metadata\
	data/escobar/escobar.taxonomy : code/escobar.batch code/escobar.R
	bash code/escobar.batch


# Get Zupancic study data and process through mothur
data/zupancic/zupancic.braycurtis.0.03.lt.ave.dist\
	data/zupancic/zupancic.groups.ave-std.summary\
	data/zupancic/zupancic.0.03.subsample.shared\
	data/zupancic/zupancic.metadata\
	data/zupancic/zupancic.taxonomy : code/zupancic.batch code/zupancic.R
	bash code/zupancic.batch


# Get Wu study data and process through mothur
data/wu/wu.braycurtis.0.03.lt.ave.dist\
	data/wu/wu.groups.ave-std.summary\
	data/wu/wu.0.03.subsample.shared\
	data/wu/wu.metadata\
	data/wu/wu.taxonomy : code/wu.batch code/wu.R
	bash code/wu.batch


# turnbaugh


# Get Goodrich study data and process through mothur

goodrich : code/goodrich.batch
	bash code/goodrich.batch
