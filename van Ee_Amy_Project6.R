# CSCB
# 27 April 2021

#--------------------------------- PROJECT 6 ----------------------------------#

#------------------------------------------------------------------------------#
#---------------------------------  SET UP  -----------------------------------#
#------------------------------------------------------------------------------#

# load file
source("helper.R")

################################### packages ################################### 
library(Seurat)
library(SeuratDisk)
library(SeuratData)


library(dplyr)
library(scmap)
library(SingleCellExperiment)
library(slingshot)
library(singleCellNet)
library(tidyverse)

# library(pheatmap)
# library(gam)
# library(RColorBrewer)
# library(rgl)
# library(singleCellNet)


##################################### data ##################################### 

# given datasets
refData <- LoadH5Seurat("TS_Ref_2021.h5seurat")
queryData <- LoadH5Seurat("Lu_Query_2021.h5seurat")

# ortholog conversion for cross species analysis
# oTab <- utils_loadObject("human_mouse_genes_Jul_24_2018.rda")
# aa = csRenameOrth(expQuery = queryData, expTrain = refData, orthTable = oTab)
# expQueryOrth <- aa[['expQuery']]
# expTrainOrth <- aa[['expTrain']]

oTab <- read.csv("oTab.csv")

newGeneNames <- c()
# for each human gene, find corresponding mouse gene
for (gene in rownames(queryData@assays$RNA@counts)){
  # gene in human
  if (gene %in% oTab[, 2]){
    # find index where gene matches
    geneInd <- which(oTab[,2] == gene)
    # add corresponding gene from mouse
    gene2add <- oTab[geneInd, 3]
  }
  else{
    gene2add <- gene
  }
  newGeneNames <- c(newGeneNames, gene2add)
}




write.csv(as.matrix(queryData@assays$RNA@counts), 'queryData.csv', sep = ',',
          row.names = F, col.names = T, quote = F)

expMat <- read.csv('queryData.csv', sep = ',', header = T)
rownames(expMat) = newGeneNames


counts <- queryData@assays$RNA@counts



# outside dataset to test
stTM <- utils_loadObject("sampTab_TM_053018.rda") # metadata
expTMraw <- utils_loadObject ("expTM_Raw_053018.rda") # expression matrix

extraData <- CreateSeuratObject(expTMraw, project = "SeuratProject", assay = "RNA",
                                min.cells = 0, min.features = 0, names.field = 1,
                                names.delim = "-", meta.data = stTM)

# query.h5@meta.data[["sample"]]


#------------------------------------------------------------------------------#
#-----------------------------------  PART I  ---------------------------------#
#------------------------------------------------------------------------------#


################################## cell-typing ################################# 

querySeurat <- queryData
refSeurat <- refData

# celltype <- function(refSeurat, querySeurat){
#----------------------------- process query --------------------------------

# normalize data
querySeurat <- NormalizeData(querySeurat, normalization.method = "LogNormalize", 
                         scale.factor = 10000)

# scale data
all.genes <- rownames(querySeurat)
querySeurat <- ScaleData(querySeurat, features = all.genes)

# perform PCA
querySeurat <- FindVariableFeatures(object = querySeurat)
querySeurat <- RunPCA(querySeurat, features = VariableFeatures(object = querySeurat))

# cluster
querySeurat <- FindNeighbors(querySeurat, dims = 1:20)
querySeurat <- FindClusters(querySeurat, resolution = 0.7)
qSeurat_clusters <- Idents(querySeurat)


#--------------------------- convert to SCE ----------------------------------

# convert Seurat to SingleCellExperiment for cell-typing
qSCE <- as.SingleCellExperiment(querySeurat)

# update
rowData(qSCE)$feature_symbol <- rownames(qSCE)



# convert Seurat to SingleCellExperiment for cell-typing
rSCE <- as.SingleCellExperiment(refSeurat)

# update
rowData(rSCE)$feature_symbol <- rownames(rSCE)

# perform feature selection
rSCE <- selectFeatures(rSCE, suppress_plot = FALSE)
# feature indexing
rSCE <- indexCluster(rSCE, cluster_col = "cell_type")

#------------------------------- cell-type -----------------------------------

# predict query cell identity by project qSCE onto indices of reference
# optimize threshold so the most cells are selected

# set up
recall_values <- c()
precision_values <- c()

bestThreshold <- 0 
currentPred <- scmapCluster(projection = qSCE, 
                            index_list = list(metadata(rSCE)$scmap_cluster_index), 
                            threshold = bestThreshold)

# for loop to find best clustering
for (i in seq(0.1, 1, 0.01)){
  newPred <- scmapCluster(projection = qSCE, 
                          index_list = list(metadata(rSCE)$scmap_cluster_index), 
                          threshold = i)
  
  if (sum(newPred$scmap_cluster_labs == "unassigned") <= 
      sum(currentPred$scmap_cluster_labs == "unassigned")){
    bestThreshold <- i
    currentPred <- newPred
  }
  
  # predictions similarities as threshold
  pred <- Val_pred[row, ]
  pred[pred >= theta] <- 1
  pred[pred != 1] <- 0
  
  # calculations for PR curve
  actual <- c(length(querySeurat$cell))
  actual[currentPred$scmap_cluster_labs == querySeurat$newAnn] <- 1
  actual[is.na(actual)] <- 0
  
  # compute values for PR
  totalActualPos <- sum(actual)
  totalActualNeg <- length(actual) - totalActualPos
  
  truePosSamples <- actual[actual == 1 & pred == 1]
  truePos <- length(truePosSamples)
  
  trueNegSamples <- actual[actual == 0 & pred == 0]
  trueNeg <- length(trueNegSamples)
  
  falseNeg <- (length(pred) - sum(pred)) - trueNeg
  falsePos <- sum(pred) - truePos
  
  # get recall
  recall <- truePos / (truePos + falseNeg)
  recall_values <- c(recall_values, recall)
  
  # get precision
  precision <- truePos / (truePos + falsePos)
  precision_values <- c(precision_values, precision)
  
}


#-------------------------- annotate each cluster -----------------------------#
querySeurat <- AddMetaData(querySeurat,
                      extractClassLabel(currentPred, qSCE),
                      col.name = 'MostFreq')


#-------------------------- plot PR curve -----------------------------#
# find AUC
AUC <- sum(diff(recall_values)*(precision_values[-1]+precision_values[-length(precision_values)]))/2


# plot final graph
plot(recall_values, precision_values, 
     width = 15, height = 6, 
     main = c("PR Curve, AUC: ", AUC), 
     xlab = "Recall", ylab = "Precision")


# return final object
# return (currentPred)
    
#


############################ evaluate performance ############################## 

# use independent dataset to test pipeline before apply to Lu

#---------------------------- subset data -------------------------------------#
# get all cell samples
numcells <- length(extraData@meta.data$cell)

# set training size to be 75% of each cell type group
smp_size <- floor(0.75 * 250)

# set random seed for reproducibility
set.seed(123)

train_ind <- sample(seq_len(numcells), size = smp_size)

# subset Bertie data into training and testing
trainData<- extraData[, train_ind]
testData <- extraData[, -train_ind]

# #---------------------------- prepare and train -------------------------------#
# # perform feature selection
# trainData <- selectFeatures(trainData, suppress_plot = FALSE)
# # feature indexing
# trainData <- indexCluster(trainData, cluster_col = "celltype")

# # train classifier by run using training data
# trainData <- celltype(refData, trainData)
# # test classifier by run using testing data
# testData <- celltype(refData, testData)

#--------------------------------- accuracy -----------------------------------#

# perform feature selection
trainData_clusters <- celltype(refData, trainData)

# find accuracy
accuracy <- (sum(trainData_clusters$scmap_cluster_labs == trainData$newAnn)) / length(trainData@meta.data$cell)



#------------------------------------ AUPR ------------------------------------#
  

# function to produce PR curve 
get_PR_point <- function(theta){
  
  
}


AUC <- 0

# vectors of x and y points for one PR graph
precision_values <- c()
recall_values <- c()

# change threshold to create curve
for (theta in seq(0, 1, by=0.001)){
  PR_point <- get_PR_point(theta, i,  row.names(Val_pred)[i])
  recall_values <- c(recall_values, PR_point[1])
  precision_values <- c(precision_values, PR_point[2])
}

# get auc and remove any NaN
recall_values <- na.omit(recall_values)
precision_values <- na.omit(precision_values)

# make sure lengths match up
if (length(recall_values) > length(precision_values)){  
  recall_values <- head(recall_values, (length(precision_values) - length(recall_values)))
} else if (length(recall_values) < length(precision_values)){  
  precision_values <- head(precision_values, (length(recall_values) - length(precision_values)))
}


# find AUC
trapezoid <- function(x,y) {
  sum(diff(x)*(y[-1]+y[-length(y)]))/2
}
AUC <- trapezoid(recall_values, precision_values) *-1


# plot final graph
plot(recall_values, precision_values, 
     width = 15, height = 6, 
     main = c("PR Curve for ", row.names(Val_pred)[i], "AUC: ", AUC), 
     xlab = "Recall", ylab = "Precision")



#------------------------------------------------------------------------------#
#----------------------------------  PART II  ---------------------------------#
#------------------------------------------------------------------------------#

################################# cell-type Lu #################################

# apply pipeline to cluster and celltype Lu data
queryData <- celltype(refData, queryData)
                         

############################### change in time #################################

# identify change in identity over time - slingshot

query_matrix <- queryData@reductions$pca@cell.embeddings[, 1:3]

lineage <- getLineages(query_matrix, queryData@meta.data$seurat_clusters, 
                       start.clus = '0', end.clus = c('7'))

curve <- getCurves(lineage)

psTime <- slingPseudotime(curve)

expDat = as.matrix(spangler[["RNA"]]@data)


# only want top 50 significantly correlated TFs
get50TFs <- function(num){
  # get pseudotimes for curve num
  pt <- psTime[which(!is.na(psTime[, num])), ]
  t <- pt[, num]
  
  # get p-values for each gene
  gpC <- gamFit(expDat, mmTFs, t)
  
  # determine 50 highest correlated genes
  top50 <- tail(sort(gpC), 50)
  
  return(top50)
  
}


# get top 50 TFs for each time curve
top50_1 <- get50TFs(1)
top50_2 <- get50TFs(2)


# heatmap 50 tfs (hm_ti)
hm_ti(expDat, names(top50_1), grps)
hm_ti(expDat, names(top50_2), grps)

################################################################################
# if not BBB ECs, what are they?





################################################################################
# specifically look at results from pipeline for rEC






