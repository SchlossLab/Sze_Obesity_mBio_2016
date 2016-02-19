# Looking for a Signal in the Noise: Revisiting Obesity and the Bacterial Microbiome

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
**Methods:** A total of 196 data sets were initially identified for inclusion in this meta-analysis.  A total of 9 data sets were included once all the search criteria were completed (controls from the GLNE07 cohort, Ross et al. 2015, Goodrich et al. 2014, Escobar et al. 2014, Zupancic et al. 2012, Human Microbiome Project, Wu et al. 2011, MetaHit, and Turnbaugh et al. 2009.  Three analysis piplines were pursued for the data.  First, each data set was analyzed separately for relationships with BMI, similar to previous reports on the re-analysis of the bacterial microbiome and obesity.  Second, a classical apporach for pooling data in a meta analysis was used to invesitigate the over all relative risk (RR) ratio from all studies of being obese based on the bacterial microbiome.  Finally, data normalization was used to investigate all data sets together for relationships between the microbiome and obesity.  
**Results:** For the first part of the analysis there was no significant difference (P > 0.05) between non-obese and obese individuals in Bacteroidetes, Firmicutes, or B/F ratio.  Shannon diversity, OTU richness, relative risk, or the NMDS of the Bray-Curtis distance matrix analyzed by PERMANOVA was significant in at least one study (P<0.05).  For the second part of the analysis there was a signifcant pooled difference in RR of obesity based on Shannon diversity (RR = 1.19, CI = 1.04, 1.36, (P-value = 0.0139)).  There was also a significant difference in RR based in the B/F ratio for the pooled data (RR = 0.86, CI = 0.75, 0.99, (P-value = 0.0296)).  
**Conclusion:** Our meta-analysis supports the  hypothesis that there is a difference between non-obese and obese individuals in Shannon diversity and B/F ratio, based on the currently avialable bacterial microbiome data included in this meta-analysis.  





*******
##### **Introduction:**  
Obesity is a growing worldwide health concern with approximately 20% of the youth (aged 2-19) in the United States classified as either overweight or obese [@ogden_prevalence_2014].  This number increases to approximately 35% in adults (aged 20 or older) and these statisitics have seen little change since 2003 [@ogden_prevalence_2014].  Traditionally the body mass index (BMI) has been used as the traditional method of classifying indviduals as non-obese or obese [@lichtash_body_2013].  However, other methods do exist, such as body adiposity index (BAI) but do not always correlate with important molecular measures of traditionally associated with obesity [@lichtash_body_2013].  Recently, there has been a lot of interest in the bacteiral microbiome and it's potential ability to modulate obesity [@brahe_can_2016; @dror_microbiota_2016].  If this modulation potential is true then the bacterial microbiome could have a enormous impact and role to play in the future treatment of obesity and helping to stem the rising obesity epidemic.    


There has been a continuous stream of studies that report that there is a link between the bacterial microbiome and obesity.  The first research article to solidify the bacterial microbiomes' role in obesity was the landmark study by Turnbaugh et al. [-@turnbaugh_obesity-associated_2006] which provided direct evidence that by modulating the microbiome one could have an effect on the obesity phenotype of mice.  This was quickly followed up by human studies that found observable differences in the bacterial microbiome between non-obese and obese individuals [@ley_microbial_2006; @turnbaugh_core_2009].  Since the publication of these studies there have been a steady stream of publications on the bacterial microbiome with associated BMI data but without clear explicit mention of any significant correlations between the bacterial community and obesity [@ross_16s_2015; @zupancic_analysis_2012; @nam_comparative_2011; arumugam_enterotypes_2011; @goodrich_human_2014; @yatsunenko_human_2012; @wu_linking_2011; @escobar_gut_2014].  However, two recent publications [@walters_meta-analyses_2014; @finucane_taxonomic_2014] re-analyzed some of the previously publishd data and came to the shocking conclusion that the human bacterial microbiome did not correlate with obesity.  One common critique of these previous studies was that they never looked at the data pooled together and that this could result in a significant result where each individual study alone would not.          

The purpose of this study is to perform an extensive meta-analysis of the bacterial microbiome and obesity by analyzing applying a more sytematic approach than has been used by previous studies in this area.  Additionaly, this study will include both seperate and pooled analysis to assess the bacterial microbiome and obesity. By including these two components this study will build on the previous two re-analysis studies and provide a stronger assesment of the direction of effect or lack there of on the bacterial microbiome and obesity.  In particular, the goal of this study is to provide additional clarity on the role of both bacterial diversity and the Bacteroidetes/Firmicutes ratio impact on obesity.

Although there have been numerous RT-PCR and qPCR studies directed at examining the bacterial microbiome and obesity [@guo_development_2008; @schwiertz_microbiota_2010].  The main method of choice to obtain the most complete picture of the overall microbial community is through 16S sequencing of a specific hypervariable region on the 16S rRNA gene.  The sequencers used to complete this task have varied depending on research group with the most popular being the Roche 454, Illumina MiSeq, and Illumina HiSeq.  For this report we limit ourselves to analyzing studies were 16S sequencing data is avaliable along with BMI information.  We also analyzed every study with a uniform approach utilizing the mothur software package for sequence processing and clean up and R statistical software for analysis.


*******
##### **Methods:**  
*Literature Review and Study Inclusion*
<br><br>
We attempted to stick to the original PRISMA guidelines for performing a meta-analysis as published by Moher et al [-@moher_preferred_2010]. A detailed description of the literature seach and the exact formula used for the initial PubMed search can be found in the online supplement along with the PRISMA flow diagram [Figure 1].  In brief, a total of 196 studies made it through the initial identification of studies.  Screening removed a total of 184 studies leaving a total of 12 eligible studies for this meta-analysis.  Upon a full reading of these specific manuscripts 2 were exculded and the reasons can be found in the online supplement.  A total of 10 studies were ultimately deemed eligible for qualitative synthesis.  One final study was removed from the quantitative analysis since no individual in the study was found to be able to be classified in the obese category based on the classical cutoffs. A total of 3 studies were from records identified by database searching [@turnbaugh_core_2009; @zupancic_analysis_2012; @escobar_gut_2014] and the other 6 were from records identified by other sources [@ross_16s_2015; @arumugam_enterotypes_2011; @goodrich_human_2014; @wu_linking_2011; @conlan_species-level_2012].  One of these studies (Baxter) has not been published at the time of this writing.  

*Subject Demographics and Study Outcomes*
<br><br>
The overall demographics for indivdiuals used in this meta-analysis split by their respective study can be found in table 1.  A summary of the study type, study population, whether the data set is published, and the outcomes reported can be found in table 2.  There is a large variation in reported outcomes with two studies finding bacterial community associations with obesity and three studies finding some association between either Firmicutes, Bacteroidetes, or the Bacteroidetes/Firmicutes ratio and obesity.  However, a total of three manuscripts found no correlations between the bacterial microbiome and obesity.  This is in addition to the two re-analysis studies published in 2014.  The cutoffs used in this study for the BMI groups were lean (BMI less than or equal to 24), overweight (BMI was greater than 24 but less than 30), obese (BMI greater than or equal to 30).  For analysis involving obese and non-obese groups the cutoffs used were a BMI less than 30 (lean/non-obese) and a BMI of equal to or greater than 30 (obese).
  
*Sequence Analysis Pipeline*
<br><br>
All data was publically available and the sequence data was downloaded from the NCBI Sequence Read Archive, MG-RAST, or the European Nucleotide Archive.  A total of 6/9 [@turnbaugh_core_2009; @ross_16s_2015; @zupancic_analysis_2012; @wu_linking_2011; @escobar_gut_2014; @conlan_species-level_2012] studies were processed by staying as close as possible to the mothur 454 standard operating procedure, which can be found at http://www.mothur.org/wiki/454_SOP.  A total of 2/9 [@goodrich_human_2014] of the studies were processed by staying as close as possible to the mothur MiSeq standard operating procedure found at http://www.mothur.org/wiki/MiSeq_SOP.  The final data set [@arumugam_enterotypes_2011] that was used took advantage of the Metaphlan2 package [@truong_metaphlan2_2015] to generate the necessary data for this analysis.  A comparison of the similarity of Metaphlan2 to 16S sequencing using the Human Microbiome Project (HMP) data set can be found in the online supplement [Figure S1-S3]. For the studies that utilized 16S rRNA gene sequencing three studies used primers targeting the V1-V3 region [@ross_16s_2015; @escobar_gut_2014; @zupancic_analysis_2012], two studies used the V4 region [@goodrich_human_2014], and one study used the V1-V2 region [@wu_linking_2011].  Where there were two choices of data sets to use we used the V3-V5 data of the HMP [@conlan_species-level_2012] and the V2 data from the Turnbaugh study [-@turnbaugh_core_2009].  A detailed walkthrough of the entire analysis pipeline can be found at https://github.com/SchlossLab/Sze_ObesityMicrobiomeMeta_PeerJ_2016/blob/master/code/sequenceProcessing.md.

*Data Analysis*
<br><br>
All analysis was performed using R version 3.2.2 and R studio version 0.99.491.  The overall analysis was split into three general categories.  The first category involved following the same method employed by Finucane et al. and Walters et al. [-@finucane_taxonomic_2014; -@walters_meta-analyses_2014] in which each study was re-analyzed seperately for associations with BMI.  The main associations with BMI that were explored in each seperate data set were Bacteroidetes phyla, Firmicutes phyla, Bacteroidetes/Firmicutes (B/F) ratio, Shannon diversity, OTU richness, NMDS visualization of Bray-Curtis distance matrix PERMANOVA, and relative risk (RR).  A two-way T-test was performed for comparison for non-obese versus obese.  An ANOVA was used to investigate all BMI categories (if available) with a Tukey post hoc test. Third, a Pearson test for linear correlation between the variables was used.  Finally, Random forest analysis to identify the best variables for classification of obese and non-obese was performed to look if there were similarities between data sets that we may have missed.  The R package AUCRF (v1.1) was used for this.
<br><br>
The second category involved pooling all the data from all the studies together to asses whether the combined direction of all the studies increased an indivdiual's RR for being obese.  The two measurements that were investigated in detail for this part were Shannon diversity and the B/F ratio.  To generate the RR the median value of each study was taken and the number of those that were non-obese or obese in each group was recorded and a Fisher exact-test was performed.  The relative risk from the the Fisher exact-test was obtained using the epiR (v0.9-6.9) package in R.  To obtain the pooled RR and confidence intervals (CI) the metafor (v1.9-8) package available for the R statistics software was utilized.  Finally, to assess potential bias we also analyzed the two pooled analysis using a funnel plot. 
<br><br>
The third and final category involved using  Zscore normalizations to pool all the data together to investigate the correlation between BMI and Shannon diversity, Bacteroidetes, Firmicutes, and B/F ratio.  Using the Zscore values a two-way t-test was used to investigate if there was a difference between non-obese and obese individuals.  Differences in BMI categories (if data was available) was investigated using ANOVA with a tukey post hoc test.  Linear correlations between these variables and BMI were explored using the Pearson test.  Where appropriate the Bonferroni correction was applied for multiple comparisons.  In all other cases a P-value under 0.05 was considered to be significant.




*******
##### Results:

*Seperate Analysis:*

We used the categorical variable of obese or non-obese to analyze each data set seperately for significant correlations to Bacteroidetes, Firmicutes, Bacteroidetes/Firmicute (B/F) ratio, Shannon diversity, OTU richness, evenness, and Bray-Curtis distance matrix which was a similar approach taken by two previous reviews [@walters_meta-analyses_2014; @finucane_taxonomic_2014].  Using a P-value cutoff for significance as less than 0.05 we found a total of 1 significant result for Shannon diversity, 3 significant results for OTU richness, and 1 signficant result for the Bray-Curtis distance matrix.  There were no significant results for Bacteroidetes, Firmicutes, or the B/F ratio.  The full summary of the exact P-values for each respective study and variable test can be found in table 3.

When investigating how well the OTUs in each study were able to classify obese and non-obese individuals the average Out Of Bag (OOB) Area Under the Curve (AUC) was 0.7553 +/- 0.0888.  The total range of the total number of variables used for the classification ranged from 3 (Ross and Escobar) to 32 (Zupancic).  The total number of trees and node size was standardized to 1000 and 20 respectively for each data set.  Perhaps surprisingly there was very little overlap in the variables used by the random forest algorithm for the classification of obese and non-obese individuals by study.  The three most common variables was for an OTU that could be classified to the family Ruminococcaceae (7/9 studies), an OTU that could be classified to the family Lachnospiraceae (5/9 studies), and the Shannon diversity metric (4/9 studies).   


*Classical Pooled Meta-analysis:*

When the Shannon diversity RR, all studies but Wu et al had those with lower Shannon diversity then the median as being at higher risk for obesity [Figure 2].  This was confirmed with the pooled analysis showing an increased RR of 1.19 (CI 1.04, 1.36) for obesity in the lower then median Shannon diversity group (P-value = 0.0139).  For the B/F ratio RR most studies showed that a low B/F ratio tended to have less obese indivdiuals than those with a high B/F ratio [Figure 3].  A notable exception to this was the HMP data set and the Baxter data set [Figure 3].  The pooled RR was 0.86 (CI 0.75, 0.99) suggesting that a low B/F ratio does indeed have less obese individuals (P-value = 0.0296).

Using a funnel plot to asses the bias of the 9 data sets used in this meta-analysis we can observe that there does not seem to be any significant  bias [Figure 4].  Both funnel plots have almost equal number of studies scattered on either side of the estimated true value.

*ZScore Normalization Pooled Analysis:*

Using this approach we found that there was a significant difference in the normalized Shannon diversity between the obese and non-obese group (P-value = 0.00493) [Figure 5A].  Further there was also a significant difference in the normalized Shannon diversity variance when the data was split between the three BMI groups (P-value = 0.03826) [Figure 5B].  Using a Tukey post hoc test this difference was driven by the overweight versus the obese group (P-value = 0.03945).  With respect to the B/F ratio using the ZScore normalization we found that there was no significant difference between obese and non-obese indivdiuals (P-value = 0.4516) [Figure 5C].  Additionally, we did not observe a difference when the indivdiuals were seperated by BMI groups (P-value = 0.3317) [Figure 5D].


*******
##### Discussion:

This study helps to provide a little more clarity to the on going debate of whether or not there are specific bacterial microbiome variables that are influenced by obesity.  Two previous reviews [@walters_meta-analyses_2014; @finucane_taxonomic_2014] have stated that the data do not support that there is a bacterial microbiome difference between non-obese and obese individuals.  However, neither study really made an attempt to pool the existing data together to try and harness the additional power that this would give.  Here we perform an extensive literature review of the existing studies on the bacterial microbiome and obesity and perform a meta-analysis on the studies that remained based on our inclusion and exclusion criteria.  We analyze the data using three different approaches.  The first approach uses the method taken by the previous reviews [@walters_meta-analyses_2014; @finucane_taxonomic_2014], the second approach takes a meta-analysis approach using RR representation by forest plots and funnel plots, the final approach also takes a meta-analysis approach but uses ZScore normalizations and does not use RR.

Overall, the first approach in this study by and large agrees with the prevous two reviews[@walters_meta-analyses_2014; @finucane_taxonomic_2014], in that for the majority of studies, no significant correlation between the bacterial microbiome and obestiy could be observed.  However, the other two approaches using the pooled data was able to identify a correlation between measures of the bacterial microbiome and obesity [Figure 2 & 3].  For Shannon diversity both approaches using the pooled data found a correlation while only one of the approaches for the B/F ratio found a correlation.  Taken together this analysis provides a preliminary finding that there is a small but detectable difference between non-obese and obese individuals with respect to both Shannon diversity and the B/F ratio.  However, our analysis would indicate that the signal for Shannon diversity may be more robust than that of the B/F ratio since both methods of pooled analysis found Shannon diversity significant while only one did for the B/F ratio.

This study does have a few limitations that are worth mentioning.  First,  it has been documented that in a number of groups BMI is not a good classification of obesity [@who_expert_consultation_appropriate_2004; @rahman_accuracy_2010].  Of particular note, the asian population has a seperate set of guidelines for obesity [@who_expert_consultation_appropriate_2004].  Although the asain population in the studies analyzed is very small [Table 1] it is possible that this may contribute to some of the noise in the data.  However, for the measures, specifically the RR, there does not seem to be a study bias associated with this [Figure 4].  We are also limited in our BMI classifications by what the previous studies published and performed.  In general, not enough studies have been done investigating the bacterial microbiome in the context of other markers for obesity or for different ethnicities to be able to perform the analysis that we did.  Thus as more specific studies are published that address ethnicity and better markers of obesity it is possible that the noise associated with these measurements will decrease.         

Second, although we make attempts to normalize the data for the pooled analysis so that we can compare the bacterial microbiome across studies it is possible that the different variable regions sampled, sequencing machines used, and the quality of the output data could have an impact on the overall conclusions reached.  However, what our analysis shows is that despite these differences, we are still able to identify an overall signal from the data.  

Third, this meta-analysis and review of the data to date is still rather preliminary.  Only a total of 1542 individuals were included in the pooled analysis which is relatively small compared to other meta-analysis in non-microbiome fields of study.  It will be necessary to continuously update this analysis, as more studies become available, in order to obtain a more accurate estimate of the true effect that obesity has on the bacterial microbiome within humans.   

Fourth, although a significant result exists for the RR for both Shannon diversity and the B/F ratio, it is small.  This is also the case for the ZScore normalization analysis of the Shannon diversity.  As an example, a lower than median Shannon diversity only has an overall 19% increase of obesity in it versus the higher than median Shannon diversity group.  So even though there is a detectable signal it remains to be seen if this is a biologically significant phenomenon that could have a tangible effect on obesity. 

Despite the many problems that could add more noise to the measurements used, this study provides a reasonable first pass attempt to pool this data together in the effort to find a correlation between the bacterial microbiome and obesity.  In doing so it provides a tanatalizing result that despite the large differences between studies it is still possible to detect a correlation between the bacterial microbiome and obesity.  However, it should be cautioned that there are still many aspects of this analysis that could be improved upon and it truly is only a first attempt at trying to summerize and analyze the pooled data that is available.  We hope that this study helps to add important information to the on-going dialogue of the bacterial microbiome's role in obesity.



*******
##### Acknowledgements:

The authors would like to thank Nielson Baxter for his help and suggestions during the sequence analysis portion.  We would also like to thank Shawn Whitefield for her help with making sure the manuscript adhered to the guidelines set out for a meta-analysis.  Finally, we would like to thank Rick Bushman, Hongzhe Li, and Pixu Shi in agreeing to release and helping to get a hold of the subject specific metadata from the Wu et al study.  


*******
##### Tables:

**Table 1:Summary Demographics of Individuals used in the Meta-analysis.**     

|     Study     |  Age (Mean +/- SD)  |  Sex (F&#124;M)  |  European Ancestry (%)  |  BMI (Mean +/- SD)  |     Min BMI     |     Max BMI     |
|:-------------:|:-------------------:|:----------------:|:-----------------------:|:-------------------:|:---------------:|:---------------:|
|    Baxter     |    54.3 +/- 9.9     |   111&#124;61    |          87.8           |    27.0 +/- 5.3     |      17.5       |      46.9       |
|     Ross      |    57.0 +/- 11.2    |    48&#124;15    |            0            |    31.6 +/- 5.3     |      22.1       |      47.9       |
|   Goodrich    |    61.5 +/- 8.9     |    505&#124;2    |      Not Available      |    26.3 +/- 4.9     |      16.2       |      44.8       |
|    Escobar    |    38.1 +/- 11.1    |    14&#124;16    |            0            |    27.4 +/- 4.5     |      19.5       |      37.6       |
|   Zupancic*   |    48.2 +/- 13.2    |   198&#124;112   |           100           |    29.2 +/- 5.2     |      16.7       |      51.1       |
|      HMP      |    26.1 +/- 5.0     |   127&#124;129   |           83            |    24.1 +/- 3.4     |       19        |       34        |
|      Wu       |    26.4 +/- 9.4     |    31&#124;27    |      Not Available      |    24.6 +/- 4.8     |       14        |      41.3       |
|   Arumugam    |    56.5 +/- 7.7     |    45&#124;40    |      Not Available      |    27.8 +/- 5.8     |      18.6       |      40.2       |
|  Turnbaugh**  |       21 - 32       |  Not Available   |          51.4           |    Not Available    |  Not Available  |  Not Available  |

*If numerical BMI information was not available then ranges and demographic information was taken from the authors published manuscript for this table.

**Only BMI group information was provided and only an age range was given.


**Table 2:Summary of Important Study Characteristics** 


|    Study    |   Study Type    |                              Population                              |  Published  |                                   Outcomes Reported                                   |
|:-----------:|:---------------:|:--------------------------------------------------------------------:|:-----------:|:-------------------------------------------------------------------------------------:|
|   Baxter    |  Observational  |  Adults without Cancer from Texas, Ontario, Michigan, Massachusetts  |     No      |                                  N/A (Unplublished)                                   |
|    Ross     |  Observational  |     Hispanic adults with and without type 2 diabetes from Texas      |     Yes     |                Mentions obesity, no specific correlations made though                 |
|  Goodrich   |  Observational  |           Adult twins and mothers from the United Kingdom            |     Yes     |    Christensenellaceae higher in low BMI individuals. No other correlations found.    |
|   Escobar   |  Observational  |                     Healthy adults from Columbia                     |     Yes     |                        Firmicutes less abundant at higher BMI                         |
|  Zupancic   |  Observational  |  Amish Adults with and without metabolic syndrome from Pennsylvania  |     Yes     |                Bacteroidetes/Firmicutes ratio adjusted BMI correlation                |
|     HMP     |  Observational  |                Healthy Adults from Texas and Missouri                |     Yes     |                              No significant correlations                              |
|     Wu      |  Observational  |                   Healthy Adults from Pennsylvania                   |     Yes     |                        Unweighted UniFrac Correlation with BMI                        |
|  Arumugam   |  Observational  |                     Healthy Adults from Denmark                      |     Yes     |                             No significant  correlations                              |
|  Turnbaugh  |  Observational  |                Adult twins and mothers from Missouri                 |     Yes     |  Decrease diversity, Bacteroidetes, and Actinobacteria correlated with increased BMI  |


**Table 3:Summary of P-values for Measurements of Interest for each Individual Study for Obese versus Normal** 


|    Study    |  Bacteroidetes  |  Firmicutes  |  B/F Ratio  |  Shannon Diversity  |  OTU Richness  |  Evenness  |  Bray Curtis  |
|:-----------:|:---------------:|:------------:|:-----------:|:-------------------:|:--------------:|:----------:|:-------------:|
|   Baxter    |      0.318      |    0.682     |    0.631    |        0.016        |     0.016      |   0.032    |     0.064     |
|    Ross     |      0.209      |    0.257     |    0.251    |        0.166        |     0.184      |   0.254    |     0.715     |
|  Goodrich   |      0.666      |    0.967     |    0.871    |        0.667        |     0.907      |   0.632    |     0.005     |
|   Escobar   |      0.093      |     0.36     |    0.466    |        0.864        |     0.195      |   0.476    |     0.089     |
|  Zupancic   |      0.958      |    0.265     |    0.99     |        0.203        |     0.035      |   0.537    |     0.573     |
|     HMP     |      0.352      |    0.636     |    0.569    |        0.585        |     0.845      |   0.356    |     0.811     |
|     Wu      |      0.917      |    0.577     |    0.674    |        0.749        |     0.751      |   0.556    |     0.892     |
|  Arumugam   |      0.191      |    0.293     |    0.053    |        0.217        |     0.988      |   0.187    |     0.055     |
|  Turnbaugh  |      0.577      |    0.407     |    0.222    |        0.07         |     0.024      |   0.125    |     0.095     |




*******
##### Figures:


![Figure 1](Sze_ObesityMicrobiomeMeta_PeerJ_2016_files/figure-html/PrismaFlowGraph.png)\


**Figure 1: PRISMA Flow Diagram of Total Records Searched.**  This was adapted from:  Moher D, Liberati A, Tetzlaff J, Altman DG, The PRISMA Group (2009). Preferred Reporting Items for Systematic Reviews and Meta-Analyses: The PRISMA Statement. PLoS Med 6(7): e1000097. doi:10.1371/journal.pmed1000097. 


![](Sze_ObesityMicrobiomeMeta_PeerJ_2016_files/figure-html/Figure_2-1.png)\

**Figure 2: Meta Analysis of the Relative Risk of Obesity Based on Shannon Diversity.**  Groups were divided for each study on high and low Shannon Diversity groups based on the median for that study.  The overall pooled relative risk was 1.19 for the low diversity group (P-value = 0.0139).


![](Sze_ObesityMicrobiomeMeta_PeerJ_2016_files/figure-html/Figure_3-1.png)\

**Figure 3: Meta Analysis of the Relative Risk of Obesity Based on Bacteroidetes/Firmicutes Ratio.**  Groups were divided for each study on high and low B/F ratio groups based on the median for that specific study.  The overall pooled relative risk was 0.86 for the low diversity group (P-value = 0.0296).


![](Sze_ObesityMicrobiomeMeta_PeerJ_2016_files/figure-html/Figure_4-1.png)\

**Figure 4: Funnel Plot of the Shannon Diversity and Bacteroidetes/Firmicutes Ratio Relative Risk.**  **A)** Overall there does not seem to be any bias associated with the studies selected for the Shannon Diversity analysis with studies falling on either side of the predicted value and those with smaller total n falling further away from this. **B)** For the Bacteroidetes/Firmicutes ratio there does not seem to be any associated bias.  The overall pattern is similar to that observed for the Shannon Diversity analysis.






![](Sze_ObesityMicrobiomeMeta_PeerJ_2016_files/figure-html/Figure_5-1.png)\


**Figure 5: Pooled Analysis of Shannon Diversity and B/F ratio by either Obesity Status or BMI Groups.**  **A)** ZScore Normalized Shannon Diversity and Obesity.  There was a significant difference between non-obese (n=1099) and obese groups (n=443) (P-value = 0.00493).  **B)** ZScore Normalized Shannon Diversity by BMI Group.  Total n for lean, overweight, and obese groups are 509, 365, 461 respectively.  Using an ANOVA with tukey post-hoc testing there was a significant diference between the overweight and obese group (P-value = 0.03944) but no difference between lean and obese group (P-value = 0.107).  **C)** ZScore Normalized B/F Ratio and Obesity.  There was no significant difference between the non-obese (n=1099) and obese (n=443) group based on B/F ratio (P-value = 0.4516).  **D)** ZScore Normalized B/F Ratio by BMI Group.  Total n for lean, overweight, and obese groups are 509, 365, 461 respectively.  There was no significant difference between any of the three groups based on ANOVA with Tukey post-hoc testing (P-value = 0.3317).  



*******
##### References:














