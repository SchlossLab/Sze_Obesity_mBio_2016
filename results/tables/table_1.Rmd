capwords <- function(s, strict = FALSE) {
		s <- as.character(s)
		cap <- function(s) paste(toupper(substring(s, 1, 1)),
									{s <- substring(s, 2); if(strict) tolower(s) else s},
														 sep = "", collapse = " " )
		sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

make_study_label <- function(dataset){
	study <- ifelse(dataset=="hmp", "HMP", capwords(dataset))
}

format_p <- function(p_value){
	char_p <- format(round(p_value, 3), nsmall=3L)

	if(p_value < 0.001){
		char_p <- "<0.001"
	}

	return(char_p)
}

get_study_summary <- function(study, beta = beta_summary){
	metadata_file <- paste0('data/', study, '/', study, '.metadata')
	metadata <- read.table(file=metadata_file, header=T, stringsAsFactors=FALSE, na.strings=c("NA", "NULL"))

	n_subjects <- nrow(metadata)

	per_obese <- format(round(100 * mean(metadata$obese), 1), nsmall=1L)

	bmi <- "NA"
	if(sum(is.na(metadata$bmi)) != nrow(metadata)){
		mean_bmi <- format(signif(mean(metadata$bmi, na.rm=T), 3), nsmall=1L)
		min_bmi <- format(signif(min(metadata$bmi, na.rm=T), 3), nsmall=1L)
		max_bmi <- format(signif(max(metadata$bmi, na.rm=T), 3), nsmall=1L)
		bmi <- paste0(mean_bmi, " (", min_bmi, "-", max_bmi, ")")
	}

	per_female <- format(round(100 * mean(tolower(metadata$sex) == 'f'), 1), nsmall=1L)

	age <- "NA"
	if(sum(is.na(metadata$age)) != nrow(metadata)){
		mean_age <- format(signif(mean(metadata$age, na.rm=T), 3), nsmall=1L)
		min_age <- format(signif(min(metadata$age, na.rm=T), 3), nsmall=1L)
		max_age <- format(signif(max(metadata$age, na.rm=T), 3), nsmall=1L)
		age <- paste0(mean_age, " (", min_age, "-", max_age, ")")
	}

	per_white <- format(round(100*mean(metadata$white), 1), nsmall=1L)

	list(study = capwords(study),
				n_subjects = n_subjects,
				per_obese = per_obese,
				bmi_summary = bmi,
				per_female = per_female,
				age_summary = age,
				per_white = per_white,
				beta_p = format_p(beta_summary[study,'p_value']))
}

beta_summary <- read.table("data/process/beta_tests.summary",
															header=T, row.names = 1)
datasets <- sort(rownames(beta_summary))

study_summary <- t(sapply(datasets, get_study_summary))

table1 <- kable(study_summary, row.names=FALSE, col.names=c("Study", "Subjects (N)", "Obese (%)", "Average BMI (Min-Max)", "Female (%)", "Average Age (Min-Max)", "Non-Hispanic White (%)", "AMOVA (P-value)"), align=c('l',rep('c', 7)), digits=3)
