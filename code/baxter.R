# We need to extract the data for people that have normal colons from the
# Baxter et al. study

shared <- read.table(file="data/baxter/baxter.shared", header=T)
metadata <- read.table(file="data/baxter/baxter.metadata", sep='\t', header=T)

normals <- as.character(metadata[metadata$dx == "normal", "sample"])

normal_shared <- shared[shared$Group %in% normals,]
write.table(normal_shared, file="data/baxter/baxter.normal.shared", quote=F, sep='\t', row.names=F)

normal_metadata <- metadata[shared$Group %in% normals,]
write.table(normal_metadata, file="data/baxter/baxter.normal.metadata", quote=F, sep='\t', row.names=F)

