library(Matrix)

tm.droplet.matrix <- readRDS("/Users/allisonhorenberg/Desktop/TM_facs_mat.rds")
counts <- readMM(file = '/Users/allisonhorenberg/Desktop/matrix.mtx')
counts <- as(counts, "dgCMatrix")
meta <- read.csv(file = '/Users/allisonhorenberg/Desktop/metadata_droplet.csv')
annotate <- read.csv(file = '/Users/allisonhorenberg/Desktop/annotations_droplet.csv')
barcodes <- read.csv(file = '/Users/allisonhorenberg/Desktop/barcodes.tsv', header = F)
barcodes <- t(as.vector(barcodes))
genes <- as.vector(rownames(tm.droplet.matrix))
rownames(counts) <- genes
colnames(counts) <- barcodes

outDat <- CreateSeuratObject(counts, project = "SeuratProject", assay = "RNA", min.cells = 0, 
  min.features = 0, names.field = 1, names.delim = "_", meta.data = NULL)

annotate[,1] <- sub("10X_P7_", "", annotate[,1])
annotate[,1] <- sub(".*_", "", annotate[,1])
for (ind in 1:length(barcodes)){
  barcodes[ind] <- substring(barcodes[ind], 1, nchar(barcodes[ind])-2)
}

barcodes <- as.vector(barcodes)

trueCellType <- c()
for (ind in 1:length(barcodes)){
  anno <- annotate[which(annotate[,1] == barcodes[ind]), 2]
  trueCellType <- c(trueCellType, anno[1])
}

trueCellType.df <- tibble("Barcodes" = barcodes, "Cell Type" = trueCellType)

write.csv(x = trueCellType.df, file = '/Users/allisonhorenberg/Desktop/TrueCellType_OutData.csv')
