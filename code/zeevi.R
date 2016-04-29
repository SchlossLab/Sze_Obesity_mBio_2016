source('code/utilities.R')
get_dependencies("dplyr")

reverse_complement <- function(string){
	paste(chartr("ATGC", "TACG", rev(unlist(strsplit(string, split="")))), collapse="")
}

get_groups_file <- function(index_file){
"data/zeevi/E13_new_R1_001.indices"
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
