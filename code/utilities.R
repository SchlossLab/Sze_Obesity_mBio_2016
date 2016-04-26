get_bacteroides_firmicutes <- function(dataset){

	shared_file <- paste0('data/', dataset, '/', dataset, '.0.03.subsample.shared')
	shared <- read.table(file=shared_file, header=T, row.names=2, stringsAsFactors=F)[,-c(1,2)]

	tax_file <- paste0('data/', dataset, '/', dataset, '.taxonomy')
	taxonomy <- read.table(file=tax_file, header=T, stringsAsFactors=F)


	n_seqs <- apply(shared, 1, sum)[[1]]

	b_otus <- taxonomy[grep("Bacteroidetes", taxonomy$Taxonomy),1]
	b_shared <- shared[,colnames(shared) %in% b_otus]
	b_count <- apply(b_shared, 1, sum)
	b_relabund <- b_count / n_seqs

	f_otus <- taxonomy[grep("Firmicutes", taxonomy$Taxonomy),1]
	f_shared <- shared[,colnames(shared) %in% f_otus]
	f_count <- apply(f_shared, 1, sum)
	f_relabund <- f_count / n_seqs

	bf_ratio <- b_relabund/f_relabund
	bf_ratio[!is.finite(bf_ratio)] <- 1e6

	return(cbind(b = b_relabund, f = f_relabund, bf = bf_ratio))
}


get_dependencies <- function(deps){

	for (dep in deps){
		if (dep %in% installed.packages()[,"Package"] == FALSE){
			install.packages(as.character(dep), quiet=TRUE);
		}
		library(dep, verbose=FALSE, character.only=TRUE)
	}

}
