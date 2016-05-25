# Looking for a Signal in the Noise: Revisiting Obesity and the Microbiome

Marc A Sze<sup>1</sup> and Patrick D Schloss<sup>1</sup>






###In-Depth Overview of Search Strategy:

The initial search strategy included looking for all papers that initially fit under the below NCBI PubMed advanced search criteria.  The terms included in this criteria were that the manuscript had to have "Bacterial Microbiome" and "Obesity, BMI, bmi, obesity" in their manuscript criteria, it was not published more than 10 years ago, they were not review articles, and it contained research on humans only.  The below formula when put into PubMed should recapitulate our initial search on the website.


```bash
(((((((((Bacterial Microbiome) AND (Obesity or bmi or body mass index or BMI or obesity) AND "last 10 years"[PDat] AND Humans[Mesh])) NOT review[ptyp]) AND "last 10 years"[PDat] AND Humans[Mesh])) AND "last 10 years"[PDat] AND Humans[Mesh])) AND "last 10 years"[PDat] AND Humans[Mesh])
```

This search yielded a total of 187 manuscripts.  From two previous other reviews of obesity and the bacterial microbiome along with knowledge of two other published papers that investigated obesity but were missed by the database search we obtained a total of 7 more articles.  We also had access to normal healthy individuals from an unpublished dataset.  This brought our total number of records to 196.  

From this total we browsed abstracts for mention of stool or feces examination, that did not involve children, was not a clinical trial for probiotics or other diet related treatments, did not only have participants with inflammatory bowel disease, the articles were in English, did not only use PCR, qPCR, or RT-PCR only for their analysis, and sequencing that used only clone libraries.  This ultimately excluded all but a total of 11 studies.

From this total of 12 studies the full text was reviewed for whether or not sequencing data was publicly available, BMI information (either categorical or continuous) was available in a supplement or, if it was not available, whether authors upon contact were willing to share this information or direct us to repositories that stored this specific information.  One study was excluded [@yatsunenko_human_2012] because it contained children and their sequencing of the 16S rRNA gene involved amplicons of only 100bp in length.  They also did not have obesity as part of their results in the actual published manuscript.  The second study was excluded because they did not have BMI information available and when contacted the authors never returned any correspondence [@zeevi_personalized_2015].  A third study was excluded since it did not use 16S rRNA gene sequencing for their bacterial microbiome analysis [@arumugam_enterotypes_2011]  

Once these 3 studies were excluded there was a total of 9 studies in the qualitative synthesis of the analysis.  Because we decided a prioi to use the standard definition for BMI group classification one study from this ten did not have any individuals who were obese by this criteria [@nam_comparative_2011] and was excluded from the final quantitative synthesis and analysis.  

*Inclusion Criteria:*

* Contains mention of Bacterial Microbiome and Obesity
* BMI, bmi, or obesity could be referenced instead of Obesity
* Not published more than 10 years ago
* Research on Humans only
* At least one specific result examining obesity and a bacterial microbiome measure
* Participants did not have Inflammatory Bowel Disease or Cancer
* Greater than 100bp single or dual end reads for 16S Sequencing
* DNA obtained from stool or feces

*Exclusion Criteria:*

* PCR, qPCR, metagenomic sequencing, or RT-PCR used as main analysis Tool
* TRFLP or clone sequencing used to asses the bacterial community
* Utilization of 100bp or less single end reads for sequencing
* Sequencing Data not publicly available for download
* BMI not available and authors do not return correspondence
* Samples were not stool or feces
* Study contained children
* Study was a review





###Supplemental Results:


**Table S1. Summary of P-values for Measurements of Interest for each Individual Study for Obese versus Normal**


|    Study    |  Bacteroidetes  |  Firmicutes  |  B/F Ratio  |  Shannon Diversity  |  OTU Richness  |  Evenness  |  Bray Curtis  |
|:-----------:|:---------------:|:------------:|:-----------:|:-------------------:|:--------------:|:----------:|:-------------:|
|   baxter    |      0.196      |    0.703     |    0.344    |        0.034        |     0.014      |   0.095    |     0.349     |
|   escobar   |      0.046      |    0.547     |    0.184    |        0.863        |     0.367      |   0.604    |     0.382     |
|     hmp     |      0.054      |    0.332     |    0.680    |        0.362        |     0.638      |   0.271    |     0.194     |
|    ross     |      0.223      |    0.452     |    0.285    |        0.201        |     0.182      |   0.305    |     0.888     |
|  schubert   |      0.182      |    0.177     |    0.171    |        0.261        |     0.862      |   0.108    |     0.637     |
|  turnbaugh  |      0.996      |    0.816     |    0.816    |        0.098        |     0.031      |   0.217    |     0.417     |
|     wu      |      0.937      |    0.776     |    0.875    |        0.647        |     0.704      |   0.552    |     0.314     |
|  zupancic   |      0.493      |    0.502     |    0.596    |        0.556        |     0.322      |   0.650    |     0.341     |
|  goodrich   |      0.617      |    0.816     |    0.861    |        0.000        |     0.000      |   0.006    |     0.003     |






**Figure S1: Boxplots of All Clostridiales OTUs Important for Non-Obese and Obese Classification.** A total of 6/8 studies had at least 1 OTU that was classified to Clostridiales.  However, overall there was a wide range in direction both between studies and within studies.




**Figure S2: Boxplots of All Lachnospiraceae OTUs Important for Non-Obese and Obese Classification.** A total of 6/8 studies had at least 1 OTU that was classified to Lachnospiraceae.  However, overall there was a wide range in direction both between studies and within studies.




**Figure S3: Boxplots of All Ruminococcaceae OTUs Important for Non-Obese and Obese Classification.** A total of 7/8 studies had at least 1 OTU that was classified to Ruminococcaceae  However, overall there was a wide range in direction both between studies and within studies.


**Table 2. Summary of BLAST alignment of Representative OTU for Clostridiales, Lachnospiraceae, and Ruminococcaceae classification OTUs**  




**Table S3. Top Sequence Similarity by Variable Region Based on the Representative OTU for Clostridiales, Lachnospiraceae, and Ruminococcaceae**  






**Figure S4: Funnel Plot of the Shannon Diversity and Bacteroidetes/Firmicutes Ratio Relative Risk.**  **A)** Overall there does not seem to be any bias associated with the studies selected for the Shannon Diversity analysis with studies falling on either side of the predicted value and those with smaller total n falling further away from this. **B)** For the Bacteroidetes/Firmicutes ratio there does not seem to be any associated bias. The overall pattern is similar to that observed for the Shannon Diversity analysis.




**Figure S5. Summary of Power and Sample Size Simulations for OTU Richness (S) for Non-obese versus Obese**  P-values listed in the legend represent the outcome of a wilcoxson rank sum test between non-obese and obese for each specific data set.  **A)** The dotted lines represent the actual effect size for each specific study.  **B)** The dotted lines represent the needed n to achieve an 80% power with the actual study effect size.



**Figure S6. Summary of Power and Sample Size Simulations for Evenness (J) for Non-obese versus Obese**  P-values listed in the legend represent the outcome of a wilcoxson rank sum test between non-obese and obese for each specific data set.  **A)** The dotted lines represent the actual effect size for each specific study.  **B)** The dotted lines represent the needed n to achieve an 80% power with the actual study effect size.




**Figure S7. Summary of Power and Sample Size Simulations for Bacteroidetes for Non-obese versus Obese**  P-values listed in the legend represent the outcome of a wilcoxson rank sum test between non-obese and obese for each specific data set.  **A)** The dotted lines represent the actual effect size for each specific study.  **B)** The dotted lines represent the needed n to achieve an 80% power with the actual study effect size.



**Figure S8. Summary of Power and Sample Size Simulations for Firmicutes for Non-obese versus Obese**  P-values listed in the legend represent the outcome of a wilcoxson rank sum test between non-obese and obese for each specific data set.  **A)** The dotted lines represent the actual effect size for each specific study.  **B)** The dotted lines represent the needed n to achieve an 80% power with the actual study effect size.




**Figure S9. Summary of Power and Sample Size Simulations for Relative Risk  for Obesity based on Shannon Diversity**  P-values listed in the legend represent the outcome of a wilcoxson rank sum test between non-obese and obese for each specific data set.  **A)** The dotted lines represent the actual effect size for each specific study.  **B)** The dotted lines represent the needed n to achieve an 80% power with the actual study effect size.



**Figure S10. Summary of Power and Sample Size Simulations for Relative Risk  for Obesity based on B/F Ratio**  P-values listed in the legend represent the outcome of a wilcoxson rank sum test between non-obese and obese for each specific data set.  **A)** The dotted lines represent the actual effect size for each specific study.  **B)** The dotted lines represent the needed n to achieve an 80% power with the actual study effect size.


*******
##### References:
