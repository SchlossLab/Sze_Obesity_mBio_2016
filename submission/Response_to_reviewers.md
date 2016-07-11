**Reviewer #1 (Comments for the Author):

Major Comments
• P.5 (line 81): The bioinformatics behind the statistical analyses are important. Include a brief summary in the main text. Specifically, it is important to clarify that you focused on studies that collected targeted 16S sequence data (not shotgun metagenomes). Also briefly mention that you did OTU clustering from scratch (per study and pooled?). Here or in the methods, clarify how alpha and beta diversity statistics were estimated. **


**• P.7 (lines 134-136): To overcome challenges with non-overlapping 16S regions, might it be possible to predict taxonomy of OTUs (which it appears you have already done for phylum level analyses) and then do a pooled community structure analysis on the taxa rather than OTUs? Alternatively, perhaps tools designed for identifying OTUs from non-overlapping shotgun marker gene reads might work. **


**• P.8 (line 162): What variability was used in the calculations (presumably that observed in the studies)? For designing future studies, it would be great to investigate the sensitivity of the results to the amount of between-sample variance, since variance does differ between studies (e.g., is higher in rural populations and is lower between moms and their newborns than between unrelated individuals). **


**• P.10 (line 204): When discussing biological relevance, it would be useful to also mention that relative abundance may not be the relevant factor for the host. Imagine an OTU (A) that produces a metabolite with a physiological effect on obesity (or whatever host trait one is studying). Suppose the number of cells and the amount of metabolite produced remains constant, but the number of cells for some other OTU (B) goes up by an order of magnitude. The relative abundance of A in the community has gone down, but this may be irrelevant to the molecular/cellular mechanism affecting obesity, especially if B has no relationship to A or obesity (e.g., it occupies a distinct niche in the gut). The same argument probably holds for alpha diversity metrics based on relative abundances - the host likely cares about specific components of the diversity not diversity per se. **


**• P.12 (line 244): Where can the pipeline be obtained? This info is in the methods but should probably be highlighted.**


**• P.12 (lines 265-267): Elaborate further on this point. The right biomarker may be a protein or pathway. Even if this is taxonomically restricted, it may be hard to detect with community-level statistics or even OTU relative abundances. In addition to the issue noted above regarding relative abundance, species-level OTUs can still miss functional variability due to large differences in gene content among strains with identical or nearly identical 16S sequences. **



**Minor Comments
• P.4 (lines 53-56): repeated text and typos: "mechanistic studies using animal models that were manipulated with antibiotics or colonization with varied communities were manipulated with antibiotics or underwent colonization with varied communities appears to support the association since these manipulation yielded" **


**P.6 (line 100) and elsewhere: When first referring to a study by name (e.g., Walters, Finucane, Goodrich) provide the citation. This info is in Table 1, but it would be helpful to provide in the text as well. **


**P.6 (line 103): In the random effects model does study have a random intercept or slope (or both)? **


**P.6 (line 88): What is the overlap of studies with the Walters and Finucane analyses (could move or copy this from the discussion)? **


**P.6 (line 88): What is the range of sample size per study (in Table 1, but would be helpful to give this overview statistic in the text)? **


**P.7 (lines 122-123): Clarify if ranges on relative risk are 95% confidence intervals or some other statistic. **


**P.7 (line 123): Fix plural: "were no significant difference" **


**P.7 (lines 132-133): Explain why directionality cannot be inferred. **


**P.13 (line 281): Is this algorithm the same as agglomerative hierarchical clustering with average linkage? Also, provide some more detail on how the read counts for OTUs were used to estimate the various alpha and beta diversity statistics. **


**P.14 (line 297): Typo: "each study within a study" **



**Reviewer #2 (Comments for the Author):

General comments:
Overall this paper is a useful meta analysis of obesity and the microbiome. I have some miscellaneous comments that are below. In addition I have one more important concern that I believe needs to be addressed.

More important concern:
Given that their analysis depends on them compiling data from multiple studies and reanalyzing that data, I believe it is necessary for these authors to make sure that all components of this analysis are reproducible. I commend the authors for their commitment to sharing code for all the work. However, I have a few concerns about the data. First, some of the data is from a personal web site. I would recommend that the authors here either repost that data to a more sustainable site (e.g., Figshare), get the original producers of the data to post it to a more sustainable site, or come up with some other solution. Again, I commend the authors for their efforts in openness and reproducibility, I am just concerned that some aspects of this paper will not be reproducible. **


**Additional concerns and comments:
L47 "individuals had a lower diversity than lean individuals (6)." Specify what kind of diversity **


**L54 "with antibiotics or colonization with varied communities were manipulated with antibiotics" Typo? Should this be "that were manipulated?" **


**L73 "Literature Review and Study Inclusion. We followed the Preferred Reporting Items for Systematic Reviews and Meta-Analyses (PRISMA) guidelines to identify studies to include in our meta-analysis (14)." Would be good to say something about whether this is considered an important approach in meta-analysis **


**L77 "we searched PubMed for original research studies that involved studying obesity and the..." Given the importance of covering the literature well in meta-analyses it would be good to say how Pubmed was searched.**


**L78 "We identified ten additional studies" How?**


**L80 "and obesity. We then manually curated the 196 studies to select." Numbers don't add up - 187 + 10 = 197. I looked through the Figure and I am not sure this explains all the numbers either.**


**L102 "other studies appeared to have the same trend, albeit the differences were not statistically significant." Quibble about wording here. If the result is not significant - is it a trend?**


**L110 "significant lower diversity than non-obese individuals; however, it is questionable whether the difference is biologically significant." Some justification for this statement is needed.**


**L126 "obese, it is questionable whether that risk is biologically or clinically relevant." Again, Some justification for this statement is needed.**


**L164 "each of the studies (Figures 6, S3-S8). Although there is no biological rationale for these effect sizes, they represent a range that is plausible." It would help to clarify what is meant here by "range that is plausible"?**


**L193 "datasets. This analysis demonstrated that the ability to reliably classify individuals as being obese based on the composition of their microbiome was limited" I think it would be good here to clarify that this is really the ability to classify based solely on microbiome data. It might be possible to do this much better if individuals were subclassified (e.g., by genotype, age, gender, etc)


**Some references could / should be added: lme4 (v.1.1-12) R package, vegan (v.2.3-5) R package, AUCRF (1.1) R package, pROC (1.8) R package, pwr (1.1-3) R package**


**Table 1 - would be good to add the references for the studies**


**Figure 1 - some of the boxes in the flow chart are cutoff.**
