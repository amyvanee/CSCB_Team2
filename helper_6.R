"Helper for Project 6 - Written (largely) by Travis Brady" 

#' loads an R object when you don't know the name
#' @param fname file
#'
#' @return variable
#'
#' @export
utils_loadObject<-function
(fname
 ### file name
){
  x<-load(fname);
  get(x);
}

"convertR: move between murine and human gene labels"

## function to swap between gene labels

convertR <- function(query = qryprdata, analog = anno, dir = "mouse"){ 
  
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

"diffMakR: check what genes are highly differentially expressed for a given cell type"

# function to subset based on cell type

subsetR <- function(data = refTrim, labels = cellTypes, type = 'oligodendrocyte'){
  ids = names(labels[labels == type])
  data[, ids]
}

subsetRQ <- function(data = qryTrim, labels = qryTimes, type = 'HUVEC'){
  ids = labels[labels == type]
  data[, ids]
}

# find (average) differential gene expression between reference type and all reference types

diffMakR <- function(avg = avgFin, cell = '', top = 1000, bot = 50){
  joint = c()
  avgdiff = sort(avg[,colnames(avg) == cell] - rowMeans(avg))
  top = names(tail(avgdiff, top))
  bottom = names(head(avgdiff, bot))
  joint = c(joint, list(top))
  joint = c(joint, list(bottom))
  names(joint) = c('high', 'low')
  joint
}

"comparR: find overlap between highly overexpressed genes in query cell and all references"

# compare expression profile of desired cell against averages for all types
# report list of overlaps

comparR <- function(avg = avgFin, query = qryFin, cell = 1, top = 1000, bot = 50){
  
  #  find differential expression between query cell and average across type
  joint = c()
  avgdiff = sort(query[,cell] - rowMeans(avg))
  top = names(tail(avgdiff, top))
  bottom = names(head(avgdiff, bot))
  joint = c(joint, list(top))
  joint = c(joint, list(bottom))
  names(joint) = c('high', 'low')
  
  
  # create a list  that has differential expression profile for all cell types
  BigDiff = c()
  for(name in colnames(avg)){
    differentiate = diffMakR(cell = name)
    BigDiff  = append(BigDiff, list(differentiate))
  }
  names(BigDiff) = colnames(avg)
  
  # check for number of overlaps in both overexpressed  and underexpressed genes
  
  Inter = c()
  for(name in names(BigDiff)){
    highinter = length(intersect(joint[['high']], BigDiff[[name]][['high']]))
    lowinter = length(intersect(joint[['low']], BigDiff[[name]][['low']]))
    inter = c(highinter, lowinter)
    Inter = append(Inter, list(inter))
    # print(joint[['high']])
    # print(name)
    # print(BigDiff[[name]][['high']])
  }
  names(Inter) = names(BigDiff)
  Inter
}

"quantifyR: ind the differences between a specific cell and all reference types"

# brute force subtraction of different expression values

quantifyR <- function(avg = avgFin, query = qryFin, cell = 1){
  
  differ = c()
  noms  = c()
  
  for(type in colnames(avg)){
    # print(type)
    dif = query[,cell] - avg[colnames(avg) == type]
    differ = append(differ, list(dif))
    noms = append(noms, type)
  }
  names(differ) = noms
  differ
}

"diffSumR: take output from quantifyR and find the sum of differences between that specific cell
and all cell types"

# raw differences aren't great

# instead  will loop through to find maximally differentially expressed genes
# use specially earmarked genes as high indicators of pluripotency
diffSumR <- function(quant){
  diffSum = c()
  for(i in 1:length(quant)){
    totdiff = sum(quant[[i]])
    diffSum = append(diffSum, totdiff)
  }
  names(diffSum) = colnames(avgFin)
  # names(which.min(diffSum))
  diffSum
}

totalR  <- function(comparison){
  Tot = c()
  for(name in names(comparison)){
    tot = sum(comparison[[name]])
    Tot = append(Tot, tot)
  }
  names(Tot) = names (comparison)
  Tot
}

extractR <- function(scores = D0_exp, index = 1, type = 'ESC'){
  typscr = scores[[index]][[type]]
  totscr = sum(scores[[index]])
  typscr/totscr
}

averageR <- function(scor = D0_exp, typ = 'ESC'){
  averages = c()
  for(i in 1:length(scor)){
    avg = extractR(scores = scor, index = i, type = typ)
    averages = append(averages, avg)
  }
  mean(averages)
}
