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
library(yardstick)

# library(pheatmap)
# library(gam)
# library(RColorBrewer)
# library(rgl)
# library(singleCellNet)


##################################### data ##################################### 

# given datasets
refData <- LoadH5Seurat("TS_Ref_2021.h5seurat")
queryData <- LoadH5Seurat("Lu_Query_2021.h5seurat")

################################################################################

# ortholog conversion for cross species analysis

# convert query to be mouse
oTab <- read.csv("oTab.csv")


convertR <- function(query, analog, dir){ 
  
  # create "notin" operator
  '%notin%' <- Negate('%in%')
  
  Conv = c()
  # if desired direction is mouse...
  if(dir == "mouse"){
    # ignore  all genes that don't have analogs listed
    query = query[rownames(query) %in% analog$human,]
    # loop through each gene in column 1 and...
    for(gene in rownames(query)){
      # replace with mouse analog
      if(gene %in% analog$human){
        conv = analog[analog$human == gene,]$mouse
        Conv = append(Conv, conv)
      }
    }
  }
  
  # if desired direction is human
  else if(dir == "human"){
    #ignore all genes that don't have analogs listed
    query = query[rownames(query) %in% analog$mouse,]
    # loop through each gene in column 1 and...
    for(gene in rownames(query)){
      # replace with human analog
      if(gene %in% analog$mouse){
        conv = analog[analog$mouse == gene,]$human
        Conv = append(Conv, conv)
      }
    }
  }
  rownames(query) = Conv
  query
}



queryExpMat <- convertR(as.matrix(queryData@assays$RNA@counts), oTab, "mouse")
queryExpMat[, 1] <- rownames(queryExpMat)

queryData <- CreateSeuratObject(queryExpMat, project = "SeuratProject", assay = "RNA",
                                min.cells = 0, min.features = 0, names.field = 1,
                                names.delim = "-", meta.data = queryData@meta.data)


# newGeneNames <- c()
# # for each human gene, find corresponding mouse gene
# for (gene in rownames(queryData@assays$RNA@counts)){
#   # gene in human
#   if (gene %in% oTab[, 2]){
#     # find index where gene matches
#     geneInd <- which(oTab[,2] == gene)
#     # add corresponding gene from mouse
#     gene2add <- oTab[geneInd, 3]
#   }
#   else{
#     gene2add <- gene
#   }
#   newGeneNames <- c(newGeneNames, gene2add)
# }
#
# already done
# write.csv(as.matrix(queryData@assays$RNA@counts), 'queryData.csv', sep = ',',
# row.names = F, col.names = T, quote = F)

# expMat <- read.csv('queryData.csv', sep = ',', header = T)
# newGeneNames[2217] = "Cebc7"
# rownames(expMat) = newGeneNames


################################################################################
# outside dataset to test (mouse)

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

celltype <- function(refSeurat, querySeurat){
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
  
  
  bestThreshold <- 0 
  currentPred <- scmapCluster(projection = qSCE, 
                              index_list = list(metadata(rSCE)$scmap_cluster_index))

  
  # for loop to find best clustering with fewest unassigned
  for (i in seq(0.1, 1, 0.01)){
    newPred <- scmapCluster(projection = qSCE, 
                            index_list = list(metadata(rSCE)$scmap_cluster_index), 
                            threshold = i)
    
    if (sum(newPred$scmap_cluster_labs == "unassigned") <= 
        sum(currentPred$scmap_cluster_labs == "unassigned")){
      bestThreshold <- i
      currentPred <- newPred
    }
  }


  # #-------------------------- annotate each cluster -----------------------------#
  # 
  clusterLabels <- extractClassLabel(currentPred, qSCE)
  names(clusterLabels) <- levels(querySeurat)
  querySeurat <- RenameIdents(querySeurat,  clusterLabels)

  
  # #-------------------------- UMAP show clusters -----------------------------#

  DimPlot(spangler, group.by = 'Idents', reduction = "umap")
  
  
  # return SCMAP output and query Seurat labelled
  return (data.frame(currentPred, querySeurat))
  
  
}


############################ evaluate performance ############################## 

# use independent dataset to test pipeline before apply to Lu

#---------------------------- subset data -------------------------------------#
# get all cell samples
numcells <- length(extraData@meta.data$cell)

# just use 50% since otherwise too large
smp_size <- floor(0.50 * 250)

# set random seed for reproducibility
set.seed(123)

train_ind <- sample(seq_len(numcells), size = smp_size)

# subset Bertie data into training and testing
extraData<- extraData[, train_ind]



# #------------------------------ pipeline ------------------------------------#
# run pipeline
extraData_clusters <- celltype(refData, extraData)[1]




#--------------------------------- accuracy -----------------------------------#

# find accuracy by compare predictions to actual labels
accuracy <- (sum(extraData_clusters$scmap_cluster_labs == extraData$newAnn)) / length(extraData@meta.data$cell)





#------------------------------------ AUPR ------------------------------------#

# PR curve
# similarities to endothelial 

# function to produce PR curve for one cellType
get_PR_point <- function(theta, similarities, trueTypes){
  
  cellType <- 'endothelial'
  
  # pred is similarities to endothelial
  pred <- similarities
  pred[pred >= theta] <- 1
  pred[pred != 1] <- 0
  
  # convert cellType classifications into binary for each sample
  actual <- c(length(trueTypes))
  actual[trueTypes == cellType] <- 1
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
  
  # get precision
  precision <- truePos / (truePos + falsePos)
  
  return(c(recall, precision))
}


# vectors of x and y points for one PR graph
precision_values <- c()
recall_values <- c()

# change threshold to create curve
for (theta in seq(0, 1, by=0.001)){
  
  # make sure input similarities to EC specifically!!
  PR_point <- get_PR_point(theta,  extraData_clusters$scmap_cluster_labs, extraData$newAnn)
  
  # get coordinate points
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
     main = c("PR Curve for Similarity to Endothelial Cell, AUPR: ", AUC, "Accuracy:", accuracy), 
     xlab = "Recall", ylab = "Precision")



#------------------------------------------------------------------------------#
#----------------------------------  PART II  ---------------------------------#
#------------------------------------------------------------------------------#

################################# cell-type Lu #################################

# apply pipeline to cluster and celltype Lu data
queryData <- celltype(refData, queryData)[2]
                         

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


# subset Lu data and apply cell typing pipeline 
queryData_rEC <- queryData[which(queryData@meta.data[["sample"]] == "IMR90_rEC")]

querydata_reEC_clusters <- celltype(refData, queryData_rEC)[1]


