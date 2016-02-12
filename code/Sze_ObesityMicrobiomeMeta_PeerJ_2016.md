# A Meta Analysis of the Bacterial Microbiome and Human Obesity

Marc A Sze<sup>1</sup> and Patrick D Schloss<sup>1</sup>


<sup>1</sup> Department of Microbiology and Immunology, University of Michigan, Ann Arbor, MI, USA  
<br>
<br>
<br>
<br>
<br>
<br>
Corresponding author: Patrick D Schloss  
*Current Address:*  
1526 MSRB I  
Ann Arbor, MI  
USA 48109-5620  
Email: pschloss@umich.edu  
Phone: 734-936-2951    
<br>
<br>
**Contributions:** MAS (planned, designed, and performed experiments, completed data analysis, wrote first draft), PDS (planned and designed the study, data interpretation, and edited subsequent drafts).  
<br>
**Running Title:** Bacterial Microbiome and Obesity
<br>





*******
##### **Abstract:**
**Rationale:** It has been well accepted within the literature that changes to the bacterial microbiome are related to obesity.  However, recent reports from two different research groups re-examining existing data suggests that there may not be a strong association between the bacterial microbiome and obesity.  These two newer studies, however, did not look at the results of pooling all the data together.  In this meta analysis, we drastically increase the number of data sets examined and investigate the pooled data to provide clearity on the association between the bacterial microbiome and obesity.  
**Methods:** A total of 11 data sets with BMI data were analyzed and processed through the standard mothur pipeline for either 454 or MiSeq data.  The data included were from Baxter et al. 2015, Ross et al. 2015, Goodrich et al. 2014, Escobar et al. 2014, Zupancic et al. 2012, Yatsunenko et al. 2012, Human Microbiome Project, Nam et al. 2011, Wu et al. 2011, MetaHit, and Turnbaugh et al. 2009.  Three analysis piplines were pursued for the data.  First, each data set was analyzed separately for relationships with BMI, similar to previous reports on the re-analysis of the bacterial microbiome and obesity.  Second, a classical apporach for pooling data in a meta analysis was used to invesitigate the over all relative risk (RR) ratio from all studies of being obese based on the bacterial microbiome.  Finally, data normalization was used to investigate all data sets together for relationships between the microbiome and obesity.  
**Results:** For the first part of the analysis there was no significant difference (P > 0.05) between non-obese and obese individuals in Bacteroidetes, Firmicutes, B/F ratio, Shannon diversity, OTU richness, relative risk, or the NMDS of the Bray-Curtis distance matrix analyzed by PERMANOVA.  For the second part of the analysis there was no pooled difference in RR of obesity based on Shannon diversity (RR = 1.05, CI = 0.90, 1.23).  There was also no difference in RR based on the B/F ratio for the poole data (RR = 0.85, CI = 0.73, 1.00).  
**Conclusion:** Our meta-analysis supports the null hypothesis that there is no difference between non-obese and obese individuals, based on the currently avialable bacterial microbiome measures.  





*******
##### **Introduction:**  
Obesity is a growing worldwide health concern with approximately 20% of the youth (aged 2-19) in the United States classified as either overweight or obese [@ogden_prevalence_2014].  This number increases to approximately 35% in adults (aged 20 or older) and these statisitics have seen little change since 2003 [@ogden_prevalence_2014].  Traditionally the body mass index (BMI) has been used as the traditional method of classifying indviduals as non-obese or obese [@lichtash_body_2013].  However, other methods do exist, such as body adiposity index (BAI) but do not always correlate with important molecular measures of traditionally associated with obesity [@lichtash_body_2013].  Recently, there has been a lot of interest in the bacteiral microbiome and it's potential ability to modulate obesity [@brahe_can_2016; @dror_microbiota_2016].  If this modulation potential is true then the bacterial microbiome could have a enormous impact and role to play in the future treatment of obesity and helping to stem the rising obesity epidemic.    


There has been a continuous stream of studies that report that there is a link between the bacterial microbiome and obesity.  The first research article to solidify the bacterial microbiomes' role in obesity was the landmark study by Turnbaugh et al. [-@turnbaugh_obesity-associated_2006] which provided direct evidence that by modulating the microbiome one could have an effect on the obesity phenotype of mice.  This was quickly followed up by human studies that found observable differences in the bacterial microbiome between non-obese and obese individuals [@ley_microbial_2006; @turnbaugh_core_2009].  Since the publication of these studies there have been a steady stream of publications on the bacterial microbiome with associated BMI data but without clear explicit mention of any significant correlations between the bacterial community and obesity [@ross_16s_2015; @zupancic_analysis_2012; @nam_comparative_2011; arumugam_enterotypes_2011; @goodrich_human_2014; @yatsunenko_human_2012; @wu_linking_2011; @escobar_gut_2014].  However, two recent publications [@walters_meta-analyses_2014; @finucane_taxonomic_2014] re-analyzed some of the previously publishd data and came to the shocking conclusion that the human bacterial microbiome did not correlate with obesity.  One common critique of these previous studies was that they never looked at the data pooled together and that this could result in a significant result where each individual study alone would not.          

The purpose of this study is to perform an extensive meta-analysis and systematic review of the bacterial microbiome and obesity by analyzing more data sets than the previous studies and including both seperate and pooled analysis to assess the bacterial microbiome and obesity. By including these two components this study will build on the previous two re-analysis studies and provide a stronger assesment of the direction of effect or lack there of on the bacterial microbiome and obesity.

Although there have been numerous RT-PCR and qPCR studies directed at examining the bacterial microbiome and obesity [@guo_development_2008; @schwiertz_microbiota_2010].  The main method of choice to obtain the most complete picture of the overall microbial community is through 16S sequencing of a specific hypervariable region on the 16S rRNA gene.  The sequencers used to complete this task have varied depending on research group with the most popular being the Roche 454, Illumina MiSeq, and Illumina HiSeq.  For this report we limit ourselves to analyzing studies were 16S sequencing data is avaliable along with BMI information.  We also analyzed every study with a uniform approach utilizing the mothur software package for sequence processing and clean up and R statistical software for analysis.


*******
##### **Methods:**  
*Literature Review and Study Inclusion*
<br><br>
A literature review was performed in pubmed for bacterial microbiome studies that also contained 16S rRNA gene sequencing data using 454, MiSeq, or HiSeq sequencing technologies.  A total of 11 studies (combined n = 2133) were identified as fullfiling this criteria.  One final study published late in 2015 used the Illumina NextSeq platform for sequencing [@zeevi_personalized_2015].  However, the BMI data for each individual was not included in the published study or the supplement.  An e-mail was sent to inquire but was never returned.  Thus it was decided not to include this data set in the subsequent meta analysis.
  
*Sequence Analysis Pipeline*
<br><br>
All data was publically available and the sequence data was downloaded from the NCBI Sequence Read Archive, MG-RAST, or the European Nucleotide Archive.  A total of 7/11 [@turnbaugh_core_2009; @ross_16s_2015; @zupancic_analysis_2012; @nam_comparative_2011; @wu_linking_2011; @escobar_gut_2014; @conlan_species-level_2012] studies were processed by staying as close as possible to the mothur 454 standard operating procedure, which can be found at http://www.mothur.org/wiki/454_SOP.  A total of 3/11 [@goodrich_human_2014; @yatsunenko_human_2012] of the studies were processed by staying as close as possible to the mothur MiSeq standard operating procedure found at http://www.mothur.org/wiki/MiSeq_SOP.  The final data set [@arumugam_enterotypes_2011] that was used took advantage of the Metaphlan2 package [@truong_metaphlan2_2015] to generate the necessary data for this analysis.  A comparison of the similarity of Metaphlan2 to 16S sequencing using the Human Microbiome Project (HMP) data set can be found in the online supplement [Figure S1 & S2].  For a few of the studies after sequence processing there was a large decrease in the number of samples available for analysis [Table 1].  A detailed walkthrough of the entire analysis pipeline can be found at https://github.com/SchlossLab/Sze_ObesityMicrobiomeMeta_PeerJ_2016/blob/master/code/sequenceProcessing.md.

*Data Analysis*
<br><br>
All analysis was performed using R version 3.2.2 and R studio version 0.99.491.  The overall analysis was split into three general categories.  The first category involved following the same method employed by Finucane et al. and Walters et al. [-@finucane_taxonomic_2014; -@walters_meta-analyses_2014] in which each study was re-analyzed seperately for associations with BMI.  The main associations with BMI that were explored in each seperate data set were Bacteroidetes phyla, Firmicutes phyla, Bacteroidetes/Firmicutes (B/F) ratio, Shannon diversity, OTU richness, NMDS visualization of Bray-Curtis distance matrix PERMANOVA, and relative risk (RR).  A two-way T-test was performed for comparison for non-obese versus obese.  An ANOVA was used to investigate all BMI categories (if available) with a Tukey post hoc test. Third, a Pearson test for linear correlation between the variables was used.  Finally, Random forest analysis to identify the best variables for classification of obese and non-obese was performed to look if there were similarities between data sets that we may have missed.  The R package AUCRF (v1.1) was used for this.
<br><br>
The second category involved pooling all the data from all the studies together to asses whether the combined direction of all the studies increased an indivdiual's RR for being obese.  The two measurements that were investigated in detail for this part were Shannon diversity and the B/F ratio.  To generate the RR the median value of each study was taken and the number of those that were non-obese or obese in each group was recorded and a Fisher exact-test was performed.  The relative risk from the the Fisher exact-test was obtained using the epiR (v0.9-6.9) package in R.  To obtain the pooled RR and confidence intervals (CI) the metafor (v1.9-8) package available for the R statistics software was utilized.  If the pooled analysis CI does not touch or cross 1.00 the pooled RR had a P-value < 0.05. 
<br><br>
The third and final category involved using  Zscore normalizations to pool all the data together to investigate the correlation between BMI and Shannon diversity, Bacteroidetes, Firmicutes, and B/F ratio.  Using the Zscore values a two-way t-test was used to investigate if there was a difference between non-obese and obese individuals.  Differences in BMI categories (if data was available) was investigated using ANOVA with a tukey post hoc test.  Linear correlations between these variables and BMI were explored using the Pearson test.  Where appropriate the Bonferroni correction was applied for multiple comparisons.  In all other cases a P-value under 0.05 was considered to be significant.




*******
##### Results:






*******
##### Discussion:







*******
**Table 1:Summary Information of Studies used in Meta-analysis.**     

|    Study     |  Year  |  Total n  |  n After Sequence Processing  |
|:------------:|:------:|:---------:|:-----------------------------:|
|    Baxter    |  2016  |    172    |              172              |
|     Ross     |  2015  |    63     |              63               |
|   Goodrich   |  2014  |    525    |              507              |
|   Escobar    |  2014  |    30     |              30               |
|   Zupancic   |  2012  |    227    |              207              |
|  Yatsunenko  |  2012  |    531    |          In Progress          |
|     HMP      |  2011  |    256    |              256              |
|     Nam      |  2011  |    18     |              18               |
|      Wu      |  2011  |    65     |              58               |
|   Arumugam   |  2010  |    85     |              85               |
|  Turnbaugh   |  2009  |    154    |              146              |




*******
##### References:














