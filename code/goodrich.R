# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

shared <- read.table("data/goodrich/goodrich.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T)

metadata_orig <- read.csv("data/goodrich/TwinsUKStudy2.csv", stringsAsFactors=F)

stopifnot(shared$Group == metadata_orig$Run_s)

sample <- metadata_orig$Run_s
sex <- ifelse(metadata_orig$sex_s == 47, 'm', 'f')
bmi <- metadata_orig$body_mass_index_s
white <- NA
age <- metadata_orig$age_s
obese <- bmi >= 30
metadata <- cbind(sample=sample, sex=sex, bmi=bmi, age=age, white=white, obese=obese)

write.table(shared, file="data/goodrich/goodrich.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/goodrich/goodrich.metadata", quote=F, sep='\t', row.names=F)
