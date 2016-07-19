STUDIES = baxter escobar hmp ross schubert turnbaugh wu zupancic goodrich zeevi
REFS = data/references
PROC = data/process
FIGS = results/figures
TABLES = results/tables

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

$(PROC)/alpha_tests.summary $(PROC)/alpha.data $(PROC)/alpha_composite.summary:\
			code/run_alpha_diversity.R code/utilities.R\
 			$(ALPHA) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

$(PROC)/beta_tests.summary : code/run_beta_diversity.R\
 			$(BETA) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

$(PROC)/relative_risk.% : code/run_relative_risk.R code/utilities.R\
 			$(ALPHA) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

$(PROC)/random_forest.otu.% : code/run_aucrf_otus.R code/cross_validate.R\
															code/utilities.R\
 			$(SHARED) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)')"

$(PROC)/random_forest.phylum.% : code/run_aucrf_taxa.R code/cross_validate.R\
																code/utilities.R\
 			$(SHARED) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)', 2)"

$(PROC)/random_forest.class.% : code/run_aucrf_taxa.R code/cross_validate.R\
																code/utilities.R\
 			$(SHARED) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)', 3)"

$(PROC)/random_forest.order.% : code/run_aucrf_taxa.R code/cross_validate.R\
																code/utilities.R\
 			$(SHARED) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)', 4)"

$(PROC)/random_forest.family.% : code/run_aucrf_taxa.R code/cross_validate.R\
																code/utilities.R\
 			$(SHARED) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)', 5)"

$(PROC)/random_forest.genus.% : code/run_aucrf_taxa.R code/cross_validate.R\
																code/utilities.R\
 			$(SHARED) $(METADATA)
	R -e "source('$<'); run('$(STUDIES)', 6)"

$(PROC)/alpha_power.% : code/run_power_analysis.R code/utilities.R\
			$(PROC)/alpha_tests.summary
	R -e "source('$<'); run_alpha('$(STUDIES)')"

$(PROC)/rr_power.% : code/run_power_analysis.R code/utilities.R\
			$(PROC)/relative_risk.summary
	R -e "source('$<'); run_rr('$(STUDIES)')"



################################################################################
#
#	Part 4: Generate figures and table
#
################################################################################

$(FIGS)/shannon_bf_ratio.tiff : code/plot_alpha_diversity.R\
												$(PROC)/alpha_composite.summary\
												$(PROC)/alpha_tests.summary
	R -e "source('$<'); build_figure(c('shannon', 'bf_ratio'))"


$(FIGS)/shannoneven_sobs_bacteroidetes_firmicutes.tiff : \
												code/plot_alpha_diversity.R\
												$(PROC)/alpha_composite.summary\
												$(PROC)/alpha_tests.summary
	R -e "source('$<'); build_figure(c('shannoneven', 'sobs', 'bacteroidetes', 'firmicutes'), leg=c(0, 10))"

$(FIGS)/funnel_plot.tiff: code/plot_funnels.R\
												data/process/relative_risk.summary
	R -e "source('$<')"

$(FIGS)/rr_shannon_bf_ratio.tiff : code/plot_rr.R\
												$(PROC)/relative_risk.composite\
												$(PROC)/relative_risk.summary
	R -e "source('$<'); build_figure(c('shannon', 'bf_ratio'))"


$(FIGS)/rr_shannoneven_sobs_bacteroidetes_firmicutes.tiff : \
												code/plot_rr.R\
												$(PROC)/relative_risk.composite\
												$(PROC)/relative_risk.summary
	R -e "source('$<'); build_figure(c('shannoneven', 'sobs', 'bacteroidetes', 'firmicutes'))"


$(FIGS)/roc_curve.tiff : code/plot_roc_data.R\
												$(PROC)/random_forest.otu.roc_data\
												$(PROC)/random_forest.otu.summary\
												$(PROC)/random_forest.genus.roc_data\
												$(PROC)/random_forest.genus.summary
	R -e "source('$<'); build_figure()"


$(FIGS)/train_test.tiff : code/plot_train_test.R\
												$(PROC)/random_forest.genus.train_test
	R -e "source('$<');"


$(FIGS)/alpha_%_power.tiff : code/plot_power.R\
												$(PROC)/alpha_power.predicted
	R -e "source('$<'); build_plots('alpha')"


$(FIGS)/rr_%_power.tiff : code/plot_power.R\
												$(PROC)/rr_power.predicted
	R -e "source('$<'); build_plots('rr')"

$(TABLES)/table_1.pdf : results/tables/table_1.Rmd $(PROC)/beta_tests.summary\
												$(METADATA)
	R -e 'render("$<")'

################################################################################
#
#	Part 5: write.paper
#
################################################################################

submission/table_1.pdf : $(TABLES)/table_1.pdf
	cp $< $@

submission/figure_1.tiff : $(FIGS)/flow_chart.png
	convert $(FIGS)/flow_chart.png submission/figure_1.tiff

submission/figure_2.tiff : $(FIGS)/shannon_bf_ratio.tiff
	cp $< $@

submission/figure_3.tiff : $(FIGS)/rr_shannon_bf_ratio.tiff
	cp $< $@

submission/figure_4.tiff : $(FIGS)/roc_curve.tiff
	cp $< $@

submission/figure_5.tiff : $(FIGS)/train_test.tiff
	cp $< $@

submission/figure_6.tiff : $(FIGS)/alpha_shannon_power.tiff
	cp $< $@


submission/supp_text.pdf : submission/supp_text.Rmd
	R -e "render('submission/supp_text.Rmd')"

submission/figure_s1.tiff : $(FIGS)/funnel_plot.tiff
	cp $< $@

submission/figure_s2.tiff : \
												$(FIGS)/shannoneven_sobs_bacteroidetes_firmicutes.tiff
	cp $< $@

submission/figure_s3.tiff : \
											$(FIGS)/rr_shannoneven_sobs_bacteroidetes_firmicutes.tiff
	cp $< $@

submission/figure_s4.tiff : $(FIGS)/alpha_bf_ratio_power.tiff
	cp $< $@

submission/figure_s5.tiff : $(FIGS)/alpha_sobs_power.tiff
	cp $< $@

submission/figure_s6.tiff : $(FIGS)/alpha_shannoneven_power.tiff
	cp $< $@

submission/figure_s7.tiff : $(FIGS)/alpha_bacteroidetes_power.tiff
	cp $< $@

submission/figure_s8.tiff : $(FIGS)/alpha_firmicutes_power.tiff
	cp $< $@

submission/figure_s9.tiff : $(FIGS)/rr_shannon_power.tiff
	cp $< $@



write.paper : submission/Sze_Obesity_mBio_2016.Rmd\
							submission/figure_1.tiff submission/figure_2.tiff\
							submission/figure_3.tiff submission/figure_4.tiff\
							submission/figure_5.tiff submission/figure_6.tiff\
							submission/table_1.pdf\
							submission/figure_s1.tiff submission/figure_s2.tiff\
							submission/figure_s3.tiff submission/figure_s4.tiff\
							submission/figure_s5.tiff submission/figure_s6.tiff\
							submission/figure_s7.tiff submission/figure_s8.tiff\
							submission/figure_s9.tiff submission/supp_text.pdf\
							$(PROC)/alpha_tests.summary $(PROC)/alpha_composite.summary\
							$(PROC)/relative_risk.summary $(PROC)/relative_risk.composite\
							$(PROC)/beta_tests.summary\
							$(PROC)/random_forest.otu.summary\
							$(PROC)/random_forest.genus.summary\
							$(PROC)/random_forest.genus.train_test\
							$(PROC)/alpha_power.predicted\
							$(PROC)/rr_power.predicted
	R -e "render('submission/Sze_Obesity_mBio_2016.Rmd', clean=FALSE)"
	mv submission/Sze_Obesity_mBio_2016.utf8.md submission/Sze_Obesity_mBio_2016.md
	rm submission/Sze_Obesity_mBio_2016.knit.md

submission/Response_to_reviewers.pdf : submission/Response_to_reviewers.md
	pandoc $< -o $@ --include-in-header=submission/header.tex

submission/Track_changes.pdf: \
					submission/Sze_Obesity_mBio_2016.md\
					submission/reference.bib\
					submission/mbio.csl\
					submission/header.tex

	OPTS="--bibliography=submission/references.bib --csl=submission/msystems.csl  --filter=pandoc-citeproc --include-in-header=submission/header.tex"
	git show 40d7145:$< > orig.md
	pandoc orig.md -o orig.tex $(OPTS)
	pandoc $< -o revised.tex $(OPTS)
	latexdiff orig.tex revised.tex > diff.tex
	pdflatex diff
	mv diff.pdf $@
	rm {revised,orig,diff}.tex
