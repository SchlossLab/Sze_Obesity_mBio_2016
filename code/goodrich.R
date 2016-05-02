# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

shared <- read.table("data/goodrich/goodrich.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T, stringsAsFactors=F)

metadata1 <- read.csv("data/goodrich/ERP006339_TwinsUK.txt", stringsAsFactors=F, sep='\t')
simple1 <- metadata1[,c("Run_s", "sex_s", "body_mass_index_s", "age_s")]


metadata2 <- read.csv("data/goodrich/ERP006342_TwinsUK.txt", stringsAsFactors=F, sep='\t')
simple2 <- metadata2[,c("Run_s", "sex_s", "body_mass_index_s", "age_s")]

simple <- rbind(simple1, simple2)

samples <- intersect(shared$Group, simple$Run_s)

shared <- shared[shared$Group %in% samples,]
shared <- shared[order(shared$Group),]

simple <- simple[simple$Run_s %in% samples,]
simple <- simple[order(simple$Run_s),]

stopifnot(shared$Group == simple$Run_s)

sample <- simple$Run_s
sex <- ifelse(simple$sex_s == 47, 'm', 'f')
bmi <- simple$body_mass_index_s
white <- NA
age <- simple$age_s
obese <- bmi >= 30
metadata <- cbind(sample=sample, sex=sex, bmi=bmi, age=age, white=white, obese=obese)

write.table(shared, file="data/goodrich/goodrich.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/goodrich/goodrich.metadata", quote=F, sep='\t', row.names=F)
