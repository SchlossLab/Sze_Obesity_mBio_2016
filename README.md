A Meta Analysis of the Bacterial Microbiome and Human Obesity
=======

##### **Abstract:**
**Rationale:** It has been well accepted within the literature that changes to the bacterial microbiome are related to obesity.  However, recent reports from two different research groups re-examining existing data suggests that there may not be a strong association between the bacterial microbiome and obesity.  These two newer studies, however, did not look at the results of pooling all the data together.  In this meta analysis, we drastically increase the number of data sets examined and investigate the pooled data to provide clearity on the association between the bacterial microbiome and obesity.  
**Methods:** A total of 196 data sets were initially identified for inclusion in this meta-analysis.  A total of 9 data sets were included once all the search criteria were completed (controls from the GLNE07 cohort, Ross et al. 2015, Goodrich et al. 2014, Escobar et al. 2014, Zupancic et al. 2012, Human Microbiome Project, Wu et al. 2011, MetaHit, and Turnbaugh et al. 2009.  Three analysis piplines were pursued for the data.  First, each data set was analyzed separately for relationships with BMI, similar to previous reports on the re-analysis of the bacterial microbiome and obesity.  Second, a classical apporach for pooling data in a meta analysis was used to invesitigate the over all relative risk (RR) ratio from all studies of being obese based on the bacterial microbiome.  Finally, data normalization was used to investigate all data sets together for relationships between the microbiome and obesity.  
**Results:** For the first part of the analysis there was no significant difference (P > 0.05) between non-obese and obese individuals in Bacteroidetes, Firmicutes, or B/F ratio.  Shannon diversity, OTU richness, relative risk, or the NMDS of the Bray-Curtis distance matrix analyzed by PERMANOVA was significant in at least one study (P<0.05).  For the second part of the analysis there was a signifcant pooled difference in RR of obesity based on Shannon diversity (RR = 1.19, CI = 1.04, 1.36, (P-value = 0.0139)).  There was also a significant difference in RR based in the B/F ratio for the pooled data (RR = 0.86, CI = 0.75, 0.99, (P-value = 0.0296)).  
**Conclusion:** Our meta-analysis supports the  hypothesis that there is a difference between non-obese and obese individuals in Shannon diversity and B/F ratio, based on the currently avialable bacterial microbiome data included in this meta-analysis.  

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

