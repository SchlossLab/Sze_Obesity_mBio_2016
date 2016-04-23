# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows of shared and metadata files must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)


shared <- read.table("data/wu/wu.trim.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T, stringsAsFactors=F)

orig_metadata <- read.table("data/wu/bmi_info.txt", header=T, stringsAsFactors=F)
orig_metadata <- orig_metadata[as.character(shared$Group), ]

sex <- c("m", "f")[orig_metadata$sex1m2f]
obese <- orig_metadata$bmi >= 30

metadata <- cbind(sample=shared$Group, sex=sex, age=orig_metadata$age, bmi=orig_metadata$bmi, obese=obese, white=NA)

write.table(metadata, file="data/wu/wu.metadata", quote=F, sep='\t', row.names=F)

write.table(shared, file="data/wu/wu.shared", quote=F, sep='\t', row.names=F)
