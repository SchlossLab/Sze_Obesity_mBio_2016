# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

all_shared <- read.table(file="data/hmp/Stool.an.shared", header=T, stringsAsFactors=F)

splits <- sapply(all_shared$Group, function(x){unlist(strsplit(x, split='\\.'))[c(1,2)]})

samples <- table(splits[1,], splits[2,])
first_sample <- apply(samples, 1, which.max)
sample_ids <- names(first_sample)
groups <- paste(names(first_sample), c("01","02","03")[first_sample], "Stool", sep='.')

shared <- all_shared[all_shared$Group %in% groups,]

stopifnot(shared$Group == groups)

write.table(shared, file="data/hmp/hmp.shared", quote=F, sep='\t', row.names=F)

categorical <- read.table(file="data/hmp/categorical.metadata", sep=' ', header=T)
continuous <- read.table(file="data/hmp/continuous.metadata", sep=' ', header=T)

stopifnot(rownames(categorical) == rownames(continuous))


good_categorical <- categorical[sample_ids,]
good_continuous <- continuous[sample_ids,]

stopifnot(rownames(categorical) == shared$Gropup)

sex <- c(Male="m", Female="f")[good_categorical$GENDER_C]
bmi <- good_continuous$DTPBMI
age <- good_continuous$AGEENR
white <- good_categorical$WHITE_C == "Yes"
obese <- bmi >= 30.0

simple_metadata <- cbind(sample=sample_ids, sex=sex, bmi=bmi, age=age, white=white, obese=obese)
colnames(simple_metadata) <- c("sample", "sex", "bmi", "age", "white", "obese")

write.table(simple_metadata, file="data/hmp/hmp.metadata", quote=F, sep='\t', row.names=F)
