source('code/utilities.R')
get_dependencies("dplyr")

reverse_complement <- function(string){
	paste(chartr("ATGC", "TACG", rev(unlist(strsplit(string, split="")))), collapse="")
}

get_groups_file <- function(index_file){
	mapping_file <- gsub("indices", "mapping", index_file)
	group_file <- gsub("indices", "groups", index_file)

	index_data <- read.table(index_file, stringsAsFactors=F)
	colnames(index_data) <- c("seq_name", "index")

	mapping_data <- read.table(mapping_file, stringsAsFactors=F, sep='\t', fill=T)
	colnames(mapping_data) <- c("sample_id",	"index", "primer", "description")

	mapping_data <- mapping_data[mapping_data$description != "",]
	mapping_data <- mapping_data[,-c(1,3)]
	mapping_data$index_rc <- sapply(mapping_data$index, reverse_complement)
	index_data$seq_name <- gsub(":", "_", index_data$seq_name)

	joined <- inner_join(index_data, mapping_data, by=c("index" = "index_rc"))
	group_data <- joined[, c("seq_name", "description")]

	write.table(group_data, group_file, row.names=F, col.names=F, quote=F, sep='\t')
}



link_shared_and_metadata <- function(){
# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

# We need to extract the data for people that have normal colons from the
# Baxter et al. study

	shared <- read.table(file="data/zeevi/zeevi.renamed.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T, stringsAsFactors=F)
	metadata <- read.table(file="data/zeevi/ENA_PRJEB11532_Metadata.csv", sep=',', header=T, stringsAsFactors=F)

	samples <- intersect(shared$Group, metadata$FD)

	shared_ordered <- shared[shared$Group %in% samples,]
	metadata_ordered <- metadata[metadata$FD %in% samples,]

	stopifnot(shared_ordered$Group == metadata_ordered$FD)



	sex <- ifelse(metadata_ordered$Gender == 0, "f", "m")
	age <- metadata_ordered$Age
	bmi <- metadata_ordered$BMI
	obese <- metadata_ordered$BMI >= 30

	metadata_simple <- data.frame(sample=metadata_ordered$FD, sex=sex, age=age, bmi=bmi, obese=obese, white=NA)

	stopifnot(shared_ordered$Group == metadata_simple$sample)

	write.table(metadata_simple, file="data/zeevi/zeevi.metadata", quote=F, sep='\t', row.names=F)

	write.table(shared_ordered, file="data/zeevi/zeevi.shared", quote=F, sep='\t', row.names=F)

}
