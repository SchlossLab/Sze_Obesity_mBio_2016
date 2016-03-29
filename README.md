A Meta Analysis of the Bacterial Microbiome and Human Obesity
=======

##### **Abstract:**
**Rationale:** The hypothesis that changes in the microbiome are associated with obesity has been widely accepted. However, two recent studies re-analyzed published data and found that when the datasets were analyzed independently they could not support the hypothesis. In this meta analysis, we retested the hypothesis by increasing the number of data sets examined and investigating the effect of pooling the results across the individual datasets.    
**Methods:** To identify datasets that we included in the meta analysis we followed the Preferred Reporting Items for Systematic Reviews and Meta-Analyses (PRISMA) guidelines we winnowed an initial set of 196 studies that tested the microbiome-obesity association hypothesis down to 8 studies that fit our criteria. We applied three analysis pipelines to the data from each study. First, we tested the association between community structure and diversity and the Bacteroides to Firmicutes ratio (B/F) with the subjects' body mass index (BMI) treated as a continuous or categorical value. We also used the random forest machine learning algorithm to identify community features that would allow us to to predict whether a person was obese. Second, we pooled the relative risk (RR) ratio of being obese based on characteristics of the microbiome for each study and tested for an association between RR and community diversity and B/F. Finally, we normalized the community diversity and B/F values within each dataset and pooled the data to test for an association with obesity.    
**Results:** When we considered each dataset in isolation, there were no significant differences (P > 0.05) in B/F between non-obese and obese individuals. Only one diversity measure (OTU richness (N=1)) was significantly different between non-obese and obese individuals (P<0.05). With individual studies, random forest was effective at identifying features within the datasets for differentiating between non-obese and obese individuals (mean AUC: 0.78; s.d.=0.08) and there was no significant difference between the studies (bonferroni adjusted P-value > 0.05). In 7/8 of the studies an OTU related to the Ruminococcaceae family was important in classification; however, the direction of the association varied between studies. When we pooled the results across studies, there was no significant difference in RR of obesity based on Shannon diversity (RR = 1.14, CI = 0.81-1.62, P-value = 0.45) and B/F ratio (RR = 1.22, CI = 0.88-1.69, P-value = 0.23). When we normalized the Shannon diversity data and pooled studies there was also no significant difference between non-obese and obese individuals (difference in mean z-score = 0.06; P = 0.33). There was also no significant difference in the normalized B/F ratio between non-obese and obese individuals (difference in mean z-score = 0.04; P = 0.44).  Finally, when completing both power and sample size calculations and simulations we found that most of the studies examined were not sufficiently powered to detect differences in any of the measurements used.   
**Conclusion:** We show that by pooling the results of many individual studies we can test the generalizability of the microbiome. Our meta-analysis provides support for the null hypothesis that there is no difference between non-obese and obese individuals in Shannon diversity and B/F ratio, that the effect sizes are small, and that a very large number of indivdiuals would have to be tested to reach sufficient study power.   

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

