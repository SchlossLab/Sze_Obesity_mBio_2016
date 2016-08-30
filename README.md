Looking for a Signal in the Noise: Revisiting Obesity and the Microbiome
=======

Two recent studies have re-analyzed published data and found that when datasets are analyzed independently there was limited support for the widely accepted hypothesis that changes in the microbiome are associated with obesity. This hypothesis was reconsidered by increasing the number of datasets and pooling the results across the individual datasets. The Preferred Reporting Items for Systematic Reviews and Meta-Analyses (PRISMA) guidelines were applied to identify 10 studies for an updated and more synthetic analysis. Alpha diversity metrics and the relative risk of obesity based on those metrics were used to identify a limited number of significant associations with obesity; however, when the results of the studies were pooled using a random effects model significant associations were observed between Shannon diversity, number of observed OTUs, and Shannon evenness and obesity status. They were not observed for the ratio of *Bacteroidetes* and *Firmicutes* or their individual relative abundances. Although these tests yielded small P-values, the difference between the Shannon diversity index of non-obese and obese individuals was 2.07%. A power analysis demonstrated that only one of the studies had sufficient power to detect a 5% difference in diversity. When Random Forest machine learning models were trained on one dataset and then tested using the other 9 datasets, the median accuracy varied between 33.01 and 64.77% (median=56.68%). Although there was support for a relationship between the microbial communities found in human feces and obesity status, this association was relatively weak and its detection is confounded by large interpersonal variation and insufficient sample sizes.

=======

### Extra Information on Data set

* To analyze full data set will need at least 200 GB of RAM
* From start to finish it will take at least 2-3 weeks to complete

=======

**Feel free to add your own data sets to the pipeline and submit a pull request to incorporate it into the full analysis**

Overview
--------

    project
    |- README          # the top level description of content
    |
    |- doc/            # documentation for the study
    |  |- notebook/    # preliminary analyses (dead branches of analysis)
    |  +- paper/       # manuscript(s), whether generated or not
    |
    |- data            # raw and primary data, are not changed once created
    |  |- references/  # reference files to be used in analysis
    |  |- raw/         # raw data, will not be altered
    |  |- mothur/      # mothur processed data
    |  +- process/     # cleaned data, will not be altered once created;
    |                  # will be committed to repo
    |
    |- code/           # any programmatic code
    |- results         # all output from workflows and analyses
    |  |- tables/      # text version of tables to be rendered with kable in R
    |  |- figures/     # graphs, likely designated for manuscript figures
    |  +- pictures/    # diagrams, images, and other non-graph graphics
    |
    |- scratch/        # temporary files that can be safely deleted or lost
    |
    |- study.Rmd       # executable Rmarkdown for this study, if applicable
    |- study.md        # Markdown (GitHub) version of the *Rmd file
    |- study.html      # HTML version of *.Rmd file
    |
    +- Makefile        # executable Makefile for this study, if applicable
