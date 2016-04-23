# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows of shared and metadata files must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)


shared <- read.table("data/turnbaugh/turnbaugh.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T, stringsAsFactors=F)

orig_metadata <- read.csv("data/turnbaugh/turnbaugh.metadata.csv", header=T, stringsAsFactors=F)

sample1_metadata <- orig_metadata[orig_metadata$Sample == 1, ]
sample1_shared <- shared[shared$Group %in% sample1_metadata$Sample_ID,]

stopifnot(sample1_shared$Group == sample1_metadata$Sample_ID)

white <- sample1_metadata$Ancestry == "EA"
obese <- sample1_metadata$BMI.category == "Obese"

metadata <- cbind(sample=sample1_shared$Group, sex="NA", age="NA", bmi="NA", obese=obese, white=white)

write.table(metadata, file="data/turnbaugh/turnbaugh.metadata", quote=F, sep='\t', row.names=F)

write.table(sample1_shared, file="data/turnbaugh/turnbaugh.shared", quote=F, sep='\t', row.names=F)
