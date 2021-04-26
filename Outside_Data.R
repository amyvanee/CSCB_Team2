library(Matrix)

tm.droplet.matrix <- readRDS("/Users/allisonhorenberg/Desktop/TM_facs_mat.rds")
counts <- readMM(file = '/Users/allisonhorenberg/Desktop/matrix.mtx')
counts <- as(counts, "dgCMatrix")
meta <- read.csv(file = '/Users/allisonhorenberg/Desktop/metadata_droplet.csv')
barcodes <- read.csv(file = '/Users/allisonhorenberg/Desktop/barcodes.tsv', header = F)
barcodes <- t(as.vector(barcodes))
genes <- as.vector(rownames(tm.droplet.matrix))
rownames(counts) <- genes
colnames(counts) <- barcodes

outDat <- CreateSeuratObject(counts, project = "SeuratProject", assay = "RNA", min.cells = 0, 
  min.features = 0, names.field = 1, names.delim = "_", meta.data = NULL)