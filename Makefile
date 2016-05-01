STUDIES = baxter escobar hmp ross schubert turnbaugh wu zupancic
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

$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align

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
#	Part 2: Get shared file, taxonomy data, representative sequence, distance
#	matrix, and metadata from individual studies
#
################################################################################

STUB = $(foreach S, $(STUDIES), data/$(S)/$(S))

ALPHA = $(addsuffix .groups.ave-std.summary,$(STUB))
BETA = $(addsuffix .braycurtis.0.03.lt.ave.dist,$(STUB))
SHARED = $(addsuffix .0.03.subsample.shared,$(STUB))
FASTA = $(addsuffix .rep.fasta,$(STUB))
TAXONOMY = $(addsuffix .taxonomy,$(STUB))
METADATA = $(addsuffix .metadata,$(STUB))

.SECONDEXPANSION:
data/%.groups.ave-std.summary\
	data/%.braycurtis.0.03.lt.ave.dist\
	data/%.0.03.subsample.shared\
	data/%.rep.fasta\
	data/%.taxonomy\
	data/%.metadata : code/$$(notdir $$*).batch code/$$(notdir $$*).R\
			$(REFS)/silva.seed.align $(REFS)/silva.v4.align\
 			$(REFS)/trainset14_032015.pds.fasta\
			$(REFS)/trainset14_032015.pds.tax
	bash $<



################################################################################
#
#	Part 3: Run alpha and beta diversity analyses, relative risk analysis, z-score
# analysis, and random forest analyses
#
################################################################################

data/process/alpha_tests.summary : code/run_alpha_diversity.R code/utilities.R\
 			$(ALPHA) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

data/process/beta_tests.summary : code/run_beta_diversity.R\
 			$(BETA) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

data/process/relative_risk.summary : code/run_relative_risk.R code/utilities.R\
 			$(ALPHA) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

data/process/z_transform.% : code/run_zscore_analysis.R code/utilities.R\
 			$(ALPHA) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

data/process/random_forest.% : code/run_aucrf_otus.R code/utilities.R\
 			$(SHARED) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

data/process/alpha_power.% : code/run_power_analysis.R code/utilities.R\
			data/process/alpha_tests.summary
	R -e "source('$<'); run_alpha()"

data/process/rr_power.% : code/run_power_analysis.R code/utilities.R\
			data/process/relative_risk.summary
	R -e "source('$<'); run_rr()"



################################################################################
#
#	Part 4: Generate figures
#
################################################################################
