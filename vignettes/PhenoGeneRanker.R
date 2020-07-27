## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PhenoGeneRankerPackage)

## -----------------------------------------------------------------------------
library(PhenoGeneRankerPackage)

## ----eval=FALSE---------------------------------------------------------------
#  walkMatrix <-CreateWalkMatrix('file.txt')
#  CreateWalkMatrix('file.txt', detectCores(), 0.4, 0.7, 0.9)

## ----eval=FALSE---------------------------------------------------------------
#  walkMatrix[[“WM”]] # accesses the walk matrix itself
#  walkMatrix[[“genes”]] # accesses the sorted gene pool nodes
#  walkMatrix[[“phenotypes”]] # accesses the sorted phenotype pool nodes
#  walkMatrix[[“gene_connectivity”]] # accesses the degree of all the genes in the network
#  walkMatrix[[“phenotype_connectivity”]] # accesses the degree of all the phenotypes in the network
#  walkMatrix[[“LG”]] # the number of genes in the network
#  walkMatrix[[“LP”]] # the number of phenotypes in the network
#  walkMatrix[[“N”]] # the number of gene pool nodes
#  walkMatrix[[“M”]] # the number of phenotype pool nodes

## ----eval=FALSE---------------------------------------------------------------
#  RWR <- RandomWalkRestart(walkMatrix, c('gene1', 'gene2'), c(), TRUE)
#  RWR <- RandomWalkRestart(CreateWalkMatrix('myFile.txt'),c('gene1'), c('phenotype1', 'phenotype2'), FALSE)
#  RWR <- RandomWalkRestart(CreateWalkMatrix('myFile.txt'),c('gene1'), c(), TRUE, 12, 0.7, 0.6, “tau”=(1,0.5,1.5), “phi”=(1,0.5,1.5))

## ----eval=FALSE---------------------------------------------------------------
#  walkMatrix <-CreateWalkMatrix('file.txt')
#  rankedGenes <-randomwalkrestart(walkMatrix, c('gene1', 'gene2'), c())
#  pValues <- RandomWalkRestart(walkMatrix ,c('gene1'), c(), TRUE, 12, 0.7, 0.6, “tau”=(1,0.5,1.5), “phi”=(1,0.5,1.5))

