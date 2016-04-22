# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

shared <- read.table("data/escobar/escobar.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T)

metadata <- read.csv("data/escobar/columbian_dataset.csv", stringsAsFactors=F)

stopifnot(shared$Group == metadata$Run_s)

sample <- metadata$Run_s
sex <- NULL
sex[grepl(" male", metadata$description_s)] <- "m"
sex[grepl("female", metadata$description_s)] <- "f"

bmi <- metadata$body_mass_index_s
white <- NA
age <- metadata$age_s
obese <- bmi >= 30
metadata <- cbind(sample=sample, sex=sex, bmi=bmi, age=age, white=white, obese=obese)

write.table(shared, file="data/escobar/escobar.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/escobar/escobar.metadata", quote=F, sep='\t', row.names=F)
