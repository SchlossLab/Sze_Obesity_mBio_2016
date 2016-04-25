library(openxlsx)

# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

# We need to extract the data for people that have normal colons from the
# schubert et al. study

shared <- read.table(file="data/schubert/clinical.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared", header=T)
metadata <- read.xlsx(xlsxFile="data/schubert/MIMARKS_cdclinical.xlsx")

normals <- as.character(metadata[metadata$disease_stat == "NonDiarrhealControl", "sample_id"])

normal_shared <- shared[shared$Group %in% normals,]
write.table(normal_shared, file="data/schubert/schubert.shared", quote=F, sep='\t', row.names=F)

normal_metadata <- metadata[metadata$sample %in% normals,]

stopifnot(normal_shared$Group == normal_metadata$sample)

#[1] "sample"       "fit_result"   "Site"         "Dx_Bin"       "dx"
#[6] "Hx_Prev"      "Hx_of_Polyps" "Age"          "Gender"       "Smoke"
#[11] "Diabetic"     "Hx_Fam_CRC"   "Height"       "Weight"       "BMI"
#[16] "White"        "Native"       "Black"        "Pacific"      "Asian"
#[21] "Other"        "Ethnic"       "NSAID"        "Abx"          "Diabetes_Med"
#[26] "stage"        "Location"

#simple_metadata <- normal_metadata[,colnames(normal_metadata) %in% keep]
sample <- normal_metadata$sample_id
white <- normal_metadata$race == "white"
sex <- tolower(normal_metadata$gender)
bmi <- NA
obese <- NA #simple_metadata$BMI >= 30
age <- normal_metadata$age
simple_metadata <- cbind(sample=sample, sex=sex, bmi=bmi, age=age, white=white, obese=obese)


write.table(simple_metadata, file="data/schubert/schubert.metadata", quote=F, sep='\t', row.names=F)
