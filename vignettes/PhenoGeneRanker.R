## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PhenoGeneRanker)

## -----------------------------------------------------------------------------
BiocManager::install("PhenoGeneRanker")
library(PhenoGeneRanker)

## ----eval=FALSE---------------------------------------------------------------
#  walkMatrix <-CreateWalkMatrix('file.txt')
#  CreateWalkMatrix('file.txt', detectCores(), 0.4, 0.7, 0.9)

## ----eval=FALSE---------------------------------------------------------------
#  walkMatrix[[“WM”]] # accesses the walk matrix itself
#  walkMatrix[[“genes”]] # sorted genes in the final network
#  walkMatrix[[“phenotypes”]] # sorted phenotypes in the final network
#  walkMatrix[[“gene_connectivity”]] # the degree of genes in the final network
#  walkMatrix[[“phenotype_connectivity”]] # the degree of phenotypes in the final network
#  walkMatrix[[“LG”]] # the number of gene layers
#  walkMatrix[[“LP”]] # the number of phenotype layers
#  walkMatrix[[“N”]] # the number of genes
#  walkMatrix[[“M”]] # the number of phenotypes

## ----eval=FALSE---------------------------------------------------------------
#  RWR <- RandomWalkRestart(walkMatrix, c('gene1', 'gene2'), c(), TRUE) # utilizes only gene seeds and generates p-values for ranks.
#  RWR <- RandomWalkRestart(CreateWalkMatrix('myFile.txt'),c('gene1'), c('phenotype1', 'phenotype2'), FALSE) # utilizes gene and phenotype seeds and does not generate p-values.
#  RWR <- RandomWalkRestart(CreateWalkMatrix('myFile.txt'),c('gene1'), c(), TRUE, 12, r=0.8, eta=0.6, tau=(1,0.5,1.5), phi=(1,0.5,1.5)) # utilizes only gene seeds, generates p-values, custom values for parameters r, eta, tau and phi.

## ----eval=FALSE---------------------------------------------------------------
#  walkMatrix <- CreateWalkMatrix('file.txt')
#  ranks <- RandomWalkRestart(walkMatrix, c('gene1', 'gene2'), c())
#  ranks_with_pvalues <- RandomWalkRestart(walkMatrix ,c('gene1', 'gene2'), c('phenotype1'), generatePValue=TRUE, numCores=12)

