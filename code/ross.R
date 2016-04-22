# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows of shared and metadata files must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)


metadata <- read.csv("data/ross/s40168-015-0072-y-s1.csv", stringsAsFactors=F)

# [1] "sampleID"          "barcode"           "proximal"
# [4] "distal"            "v_region"          "HbA1C"
# [7] "sex"               "age_at_visit"      "height"
#[10] "weight"            "BMI"               "sequence_run"
#[13] "ALT"               "HOMA"              "experiment_center"
#[16] "sample_type"       "Description"       "waist"
#[19] "hip"               "arm_circ"          "pulse2"
#[22] "labid"             "cholesterol"       "triglycerides"
#[25] "hdl_chol"          "ldl_chol"          "crp"
#[28] "mmol_gluc"         "homa_ir"           "homa_beta"

keep <- c("sampleID", "sex", "BMI", "age_at_visit")
md_keep <- metadata[,colnames(metadata) %in% keep]
md_keep$obese <- md_keep$BMI >= 34.0
md_keep$white <- FALSE

lookup <- read.csv("data/ross/Hispanic_dataset.csv", stringsAsFactors=F)
sample <- lookup[,"Run_s"]
names(sample) <- lookup$Sample_Name_s

md_keep$sampleID <- sample[md_keep$sampleID]
colnames(md_keep) <- c("sample", "sex", "age", "bmi", "obese", "white")
md_keep <- md_keep[order(md_keep$sample), ]

shared <- read.table("data/ross/ross.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T)
shared <- shared[order(shared$Group), ]

stopifnot(md_keep$sample == shared$Group)


write.table(md_keep, file="data/ross/ross.metadata", quote=F, sep='\t', row.names=F)

write.table(shared, file="data/ross/ross.shared", quote=F, sep='\t', row.names=F)
