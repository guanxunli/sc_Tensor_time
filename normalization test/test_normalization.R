library(Matrix)
source("utility.R")

#### Z normalization method
method <- 'Z-normalization'
nRun <- 10
nCells <- c(10, 50, 100, 500, 1000, 2000, 3000)
nRun <- seq_len(nRun)
nCells <- nCells

experimentDesign <- expand.grid(method, nRun, nCells)
experimentDesign <- as.data.frame.array(experimentDesign)

inputData <- readMM('data/CTL.mtx')
rownames(inputData) <- readLines('data/geneList.txt')

## generate networks
set.seed(1)
apply(experimentDesign,1,function(X){
  newFile <- paste0('networks/', method, '/', X[1], '_run', as.character(X[2]), '_',as.character(X[3]), 'cells', collapse = '')
  newFile <- gsub('[[:space:]]+','',newFile)
  for (i in 1:100){
    sCells <- sample(seq_len(ncol(inputData)), size = as.numeric(X[3]), replace = TRUE)
    tData <- inputData[,sCells]
    gList <- rownames(tData)
    nGenes <- length(gList)
    outM <- matrix(0, nrow = nGenes, ncol = nGenes)
    rownames(outM) <- colnames(outM) <- gList
    tData <- tData[apply(tData!=0,1,sum) > 0,]
    tData <- as.matrix(tData)
    if (min(colSums(tData)) > 0) break
  }
  tData <- scTenifoldNet::pcNet(tData)
  tData <- round(tData,2)
  outM[rownames(tData),colnames(tData)] <- as.matrix(tData)
  outM <- as(outM, 'dgCMatrix')
  writeMM(outM, newFile)
})

#### log normalization method
method <- 'log-normalization'
nRun <- 10
nCells <- c(10, 50, 100, 500, 1000, 2000, 3000)
nRun <- seq_len(nRun)
nCells <- nCells

experimentDesign <- expand.grid(method, nRun, nCells)
experimentDesign <- as.data.frame.array(experimentDesign)

inputData <- readMM('data/CTL.mtx')
rownames(inputData) <- readLines('data/geneList.txt')

## generate networks
set.seed(1)
apply(experimentDesign,1,function(X){
  newFile <- paste0('networks/', method, '/', X[1], '_run', as.character(X[2]), '_',as.character(X[3]), 'cells', collapse = '')
  newFile <- gsub('[[:space:]]+','',newFile)
  for (i in 1:100){
    sCells <- sample(seq_len(ncol(inputData)), size = as.numeric(X[3]), replace = TRUE)
    tData <- inputData[,sCells]
    gList <- rownames(tData)
    nGenes <- length(gList)
    outM <- matrix(0, nrow = nGenes, ncol = nGenes)
    rownames(outM) <- colnames(outM) <- gList
    tData <- tData[apply(tData!=0,1,sum) > 0,]
    tData <- as.matrix(tData)
    if (min(colSums(tData)) > 0) break
  }
  tData <- log(scTenifoldNet::cpmNormalization(tData) + 1)
  tData <- as.matrix(tData)
  tData <- pcNet(tData)
  tData <- round(tData,2)
  outM[rownames(tData),colnames(tData)] <- as.matrix(tData)
  outM <- as(outM, 'dgCMatrix')
  writeMM(outM, newFile)
})

#### new normalization method
method <- 'new-normalization'
nRun <- 10
nCells <- c(10, 50, 100, 500, 1000, 2000, 3000)
nRun <- seq_len(nRun)
nCells <- nCells

experimentDesign <- expand.grid(method, nRun, nCells)
experimentDesign <- as.data.frame.array(experimentDesign)

inputData <- readMM('data/CTL.mtx')
rownames(inputData) <- readLines('data/geneList.txt')

## generate networks
set.seed(1)
apply(experimentDesign,1,function(X){
  newFile <- paste0('networks/', method, '/', X[1], '_run', as.character(X[2]), '_',as.character(X[3]), 'cells', collapse = '')
  newFile <- gsub('[[:space:]]+','',newFile)
  for (i in 1:100){
    sCells <- sample(seq_len(ncol(inputData)), size = as.numeric(X[3]), replace = TRUE)
    tData <- inputData[,sCells]
    gList <- rownames(tData)
    nGenes <- length(gList)
    outM <- matrix(0, nrow = nGenes, ncol = nGenes)
    rownames(outM) <- colnames(outM) <- gList
    tData <- tData[apply(tData!=0,1,sum) > 0,]
    tData <- as.matrix(tData)
    if (min(colSums(tData)) > 0) break
  }
  tData <- new_Normalization(tData)
  tData <- pcNet(tData)
  tData <- round(tData,2)
  outM[rownames(tData),colnames(tData)] <- as.matrix(tData)
  outM <- as(outM, 'dgCMatrix')
  writeMM(outM, newFile)
})

#### original normalization method
method <- 'ori-normalization'
nRun <- 10
nCells <- c(10, 50, 100, 500, 1000, 2000, 3000)
nRun <- seq_len(nRun)
nCells <- nCells

experimentDesign <- expand.grid(method, nRun, nCells)
experimentDesign <- as.data.frame.array(experimentDesign)

inputData <- readMM('data/CTL.mtx')
rownames(inputData) <- readLines('data/geneList.txt')

## generate networks
set.seed(1)
apply(experimentDesign,1,function(X){
  newFile <- paste0('networks/', method, '/', X[1], '_run', as.character(X[2]), '_',as.character(X[3]), 'cells', collapse = '')
  newFile <- gsub('[[:space:]]+','',newFile)
  for (i in 1:100){
    sCells <- sample(seq_len(ncol(inputData)), size = as.numeric(X[3]), replace = TRUE)
    tData <- inputData[,sCells]
    gList <- rownames(tData)
    nGenes <- length(gList)
    outM <- matrix(0, nrow = nGenes, ncol = nGenes)
    rownames(outM) <- colnames(outM) <- gList
    tData <- tData[apply(tData!=0,1,sum) > 0,]
    tData <- as.matrix(tData)
    if (min(colSums(tData)) > 0) break
  }
  tData <- scTenifoldNet::cpmNormalization(tData)
  tData <- scTenifoldNet::pcNet(tData)
  tData <- round(tData,2)
  outM[rownames(tData),colnames(tData)] <- as.matrix(tData)
  outM <- as(outM, 'dgCMatrix')
  writeMM(outM, newFile)
})




