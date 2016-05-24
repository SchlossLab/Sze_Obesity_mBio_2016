source('code/utilities.R')
get_dependencies('vegan')

read_lt_matrix <- function(dist_file){

	file_content <- scan(dist_file, what="", quiet=T)

	n_seqs <- as.numeric(file_content[1])
	file_content <- file_content[-1]

	#insure that this is a lower triangle distance matrix
	stopifnot(n_seqs + n_seqs*(n_seqs-1)/2 == length(file_content))

	dist <- matrix(0, nrow=n_seqs, ncol=n_seqs)
	sample_names <- rep("", n_seqs)

	for(i in 1:n_seqs){
			sample_names[i] <- file_content[1]
			file_content <- file_content[-1]

			if(i > 1){
				dist[i, 1:(i-1)] <- as.numeric(file_content[1:(i-1)])
				file_content <- file_content[-(1:(i-1))]
			}
	}

	dist <- dist + t(dist)
	rownames(dist) <- sample_names

	return(dist)
}


run <- function(datasets){

	datasets <- unlist(strsplit(datasets, split=" "))

	summary_data <- NULL

	for(d in datasets){
		set.seed(1976)
		beta_file <- paste0('data/', d, '/', d, '.braycurtis.0.03.lt.ave.dist')
		beta <- read_lt_matrix(beta_file)

		metadata_file <- paste0('data/', d, '/', d, '.metadata')
		metadata <- read.table(file=metadata_file, header=T, stringsAsFactors=F)
		metadata <- metadata[metadata$sample %in% rownames(beta),]

		stopifnot(rownames(beta) == metadata$sample)

		na_obesity <- is.na(metadata$obese)
		metadata <- metadata[!na_obesity,]
		beta <- beta[!na_obesity,!na_obesity]

		stopifnot(rownames(beta) == metadata$sample)

		test <- adonis(beta~metadata$obese, permutations=9999)
		p_value <- test["aov.tab"][[1]]$'Pr(>F)'[1]

		study_data <- cbind(dataset = d, p_value = p_value)

		summary_data <- rbind(summary_data, study_data)
	}

	write.table(summary_data, file="data/process/beta_tests.summary", quote=F, sep='\t', row.names=F)

}
