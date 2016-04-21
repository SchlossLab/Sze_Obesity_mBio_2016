# Output shared file that has the same rows as the metadata file. Requirements:
#		* Rows must be in the same order
#		* Metadata file must contain sample id, gender (m/f), bmi, age,
#			white(logical), obese(logical)

# We need to extract the data for people that have normal colons from the
# Baxter et al. study

shared <- read.table(file="data/baxter/glne007.final.an.unique_list.shared", header=T)
metadata <- read.table(file="data/baxter/metadata.tsv", sep='\t', header=T)

normals <- as.character(metadata[metadata$dx == "normal", "sample"])

normal_shared <- shared[shared$Group %in% normals,]
write.table(normal_shared, file="data/baxter/baxter.shared", quote=F, sep='\t', row.names=F)

normal_metadata <- metadata[metadata$sample %in% normals,]

stopifnot(normal_shared$Group == normal_metadata$sample)

#[1] "sample"       "fit_result"   "Site"         "Dx_Bin"       "dx"
#[6] "Hx_Prev"      "Hx_of_Polyps" "Age"          "Gender"       "Smoke"
#[11] "Diabetic"     "Hx_Fam_CRC"   "Height"       "Weight"       "BMI"
#[16] "White"        "Native"       "Black"        "Pacific"      "Asian"
#[21] "Other"        "Ethnic"       "NSAID"        "Abx"          "Diabetes_Med"
#[26] "stage"        "Location"

keep <- c("sample", "Gender", "BMI", "Age", "White")

simple_metadata <- normal_metadata[,colnames(normal_metadata) %in% keep]
simple_metadata$obese <- simple_metadata$BMI >= 30
simple_metadata$White <- as.logical(simple_metadata$White)
colnames(simple_metadata) <- tolower(colnames(simple_metadata))

write.table(simple_metadata, file="data/baxter/baxter.metadata", quote=F, sep='\t', row.names=F)
