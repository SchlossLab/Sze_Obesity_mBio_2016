# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows of shared and metadata files must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

shared <- read.table("data/zupancic/zupancic.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T, stringsAsFactors=F)
#shared <- shared[order(shared$Group), ]

categorical <- read.csv("data/zupancic/amish_obesity_table2.csv", stringsAsFactors=F)
categorical <- categorical[!duplicated(categorical$submitted_sample_id_s), ]
categorical <- categorical[categorical$submitted_sample_id_s %in% shared$Group,]
categorical <- categorical[order(categorical$submitted_sample_id_s),]

continuous <- read.csv("data/zupancic/amish.metadata.csv", stringsAsFactors=F)
continuous <- continuous[!duplicated(continuous$SUBJID), ]
continuous <- continuous[continuous$SAMPID %in% shared$Group,]
continuous <- continuous[order(continuous$SAMPID),]

stopifnot(continuous$SAMPID == categorical$submitted_sample_id_s)
stopifnot(continuous$SAMPID == shared$Group)

sex <- ifelse(categorical == "female", "f", "m")
#age <- continuous$AGE
#bmi <- continuous$BMI
obese <- continuous$BMI >= 30

metadata <- cbind(sample=shared$Group, sex=sex, age=continuous$AGE, bmi=continuous$BMI, obese=obese, white=TRUE)

write.table(metadata, file="data/zupancic/zupancic.metadata", quote=F, sep='\t', row.names=F)

file.copy(from="data/zupancic/zupancic.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", to="data/zupancic/zupancic.shared", overwrite=TRUE)
