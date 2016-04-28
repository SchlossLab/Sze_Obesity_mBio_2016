library(openxlsx)

# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

# We need to extract the data for people that have normal colons from the
# schubert et al. study

shared <- read.table(file="data/schubert/clinical.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared", header=T)

metadata_orig <- read.xlsx(xlsxFile="data/schubert/MIMARKS_cdclinical.xlsx")
metadata_bmi <- read.xlsx(xlsxFile="data/schubert/ERIN_Outpatient_heights_deidentified.xlsx", rowNames=1)

normals <- as.character(metadata_orig[metadata_orig$disease_stat == "NonDiarrhealControl", "sample_id"])

w_bmi <- rownames(metadata_bmi[!is.na(metadata_bmi$BMI), ])

good_samples <- intersect(normals, w_bmi)


normal_shared <- shared[shared$Group %in% good_samples,]

write.table(normal_shared, file="data/schubert/schubert.shared", quote=F, sep='\t', row.names=F)

bmi_table <- metadata_bmi[good_samples,]
metadata_table <- metadata_orig[metadata_orig$sample %in% good_samples,]

stopifnot(normal_shared$Group == metadata_table$sample)
stopifnot(rownames(bmi_table) == metadata_table$sample)

#[1] "sample"       "fit_result"   "Site"         "Dx_Bin"       "dx"
#[6] "Hx_Prev"      "Hx_of_Polyps" "Age"          "Gender"       "Smoke"
#[11] "Diabetic"     "Hx_Fam_CRC"   "Height"       "Weight"       "BMI"
#[16] "White"        "Native"       "Black"        "Pacific"      "Asian"
#[21] "Other"        "Ethnic"       "NSAID"        "Abx"          "Diabetes_Med"
#[26] "stage"        "Location"

sample <- metadata_table$sample_id
white <- metadata_table$race == "white"
sex <- tolower(metadata_table$gender)
bmi <- bmi_table$BMI
obese <- bmi >= 30
age <- metadata_table$age
simple_metadata <- cbind(sample=sample, sex=sex, bmi=bmi, age=age, white=white, obese=obese)


write.table(simple_metadata, file="data/schubert/schubert.metadata", quote=F, sep='\t', row.names=F)
