level_regexp <- "^0\\.03"

extra_names_file <- "data/zeevi/zeevi.renamed.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table.extra.temp"
extras <- read.table(file=extra_names_file, header=T, stringsAsFactors=F)[,1]

list_files <- list.files('data/zeevi', pattern='*.list', full.names=T)
file_number <- as.numeric(gsub(".*fasta.(\\d*)\\..*", "\\1", list_files))
file_order <- order(file_number)
list_files <- list_files[file_order]

combined_data <- rep("", length(list_files))

for(i in 0:(length(list_files)-1)){

	list_line <- ""
	list_file <- scan(list_files[i+1], what=character(), sep="\n", quiet=T)
	hits <- grepl(level_regexp, list_file)

	if(sum(hits) == 1){
		list_line <- list_file[hits]
	} else {
		list_file[1] <- gsub("^unique", "0.00", list_file[1])
		cutoffs <- as.numeric(gsub("^(0\\.\\d\\d)\\t.*", "\\1", list_file))

		best_cutoff <- which.max(cutoffs[cutoffs <= 0.03])
		list_line <- list_file[best_cutoff]
	}

	combined_data[i+1] <- gsub("^0\\.\\d\\d\\t\\d*\\t", "", list_line)
}

clustered <- paste(combined_data, collapse='\t')
extra_merged <- paste(extras, collapse='\t')

combined_data <- paste(clustered, extra_merged, sep='\t')

n_otus <- sum(nchar(combined_data) - nchar(gsub("\t", "", combined_data)) + 1)

new_list_line <- paste("0.03", n_otus, combined_data, sep='\t')

write(new_list_line, file="data/zeevi/zeevi.renamed.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list")
