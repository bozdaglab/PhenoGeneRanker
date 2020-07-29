ReadSettings <- function(input.file){

  files <-  read.table(input.file, header=TRUE,sep="\t", 
                       stringsAsFactors = FALSE)
 
  LG <- sum(files$type=="gene")
  LP <- sum(files$type=="phenotype")
  output=list(
    FilesDF=files,
    LG=LG,
    LP=LP)
}
ReadNetworkLayers <- function(filesDF){
  # read the gene layer names
  filesDF <- filesDF[filesDF$type%in%c("gene", "phenotype", "bipartite"),]
  NetworkLayers <- vector("list",nrow(filesDF))
  j <- 1
  for(f in filesDF$file_name){
    NetworkLayers[[j]][["DF"]] <-  read.table(f, header=TRUE, sep="\t", 
                                              stringsAsFactors = FALSE)
    NetworkLayers[[j]][["DF"]]$from <- as.character(
                                    NetworkLayers[[j]][["DF"]]$from)
    NetworkLayers[[j]][["DF"]]$to <- as.character(NetworkLayers[[j]][["DF"]]$to)
    NetworkLayers[[j]][["type"]] <- filesDF$type[j]
    
    NetworkLayers[[j]][["graph"]]  <- graph.data.frame(
                                  NetworkLayers[[j]][["DF"]], directed = FALSE)
    
    NetworkLayers[[j]][["name"]] <- f
    j <- j+1
  }
  
  gene_pool_nodes_sorted <- GeneratePoolNodes(NetworkLayers, type = "gene")
  phenotype_pool_nodes_sorted <- GeneratePoolNodes(NetworkLayers, 
                                                   type = "phenotype")
  
  
  idx <- which(lapply(NetworkLayers, `[[`, "type") == "gene")
  for(i in idx){
    NetworkLayers[[i]][["graph"]] <- AddMissingNodesToGraph(
                                      NetworkLayers[[i]][["graph"]],
                                      "gene", gene_pool_nodes_sorted)  
  }
  
  idx <- which(lapply(NetworkLayers, `[[`, "type") == "phenotype")
  for(i in idx){
    NetworkLayers[[i]][["graph"]] <- AddMissingNodesToGraph(
                                      NetworkLayers[[i]][["graph"]],
                                      "phenotype", phenotype_pool_nodes_sorted)  
  }  
  
  output=list(
    NetworkLayers=NetworkLayers,
    gene_pool_nodes_sorted=gene_pool_nodes_sorted,
    phenotype_pool_nodes_sorted=phenotype_pool_nodes_sorted)    
  return(output)
}

ReadSeeds <- function(fileName, gene_pool_nodes_sorted, 
                      phenotype_pool_nodes_sorted){
  AllSeeds <- read.csv(fileName, header=FALSE, sep="\t", 
                       stringsAsFactors = FALSE)
  AllSeeds <- AllSeeds$V1
  SeedList <- CheckSeeds(AllSeeds,gene_pool_nodes_sorted,
                         phenotype_pool_nodes_sorted)
  return(SeedList)
}
GeneratePoolNodes <- function(FullNet, type){
  idx <- which(lapply(FullNet, `[[`, "type") == type)
  DFs <- lapply(FullNet[idx], `[[`, "DF")
  Node_Names_all <- unique(c(unlist(lapply(DFs, '[[', 'from')), 
                             unlist(lapply(DFs, '[[', 'to'))))
  ## We remove duplicates and sort
  pool_nodes_sorted <- sort(Node_Names_all)
  
  return(pool_nodes_sorted)
}
AddMissingNodesToGraph <- function(g, type, pns){
  #pns is pool_nodes_sorted
  ## We add to each layer the missing nodes of the total set of nodes, 
  #of the pool of nodes.
  Node_Names_Layer <- V(g)$name
  Missing_Nodes <- pns[which(!pns %in% Node_Names_Layer)]
  g <- add_vertices(g ,length(Missing_Nodes), name=Missing_Nodes)
  return(g)
}

AddParameters <- function(delta, zeta, lambda){
  
  if (delta > 1 || delta< 0){ 
    stop("Incorrect delta, it must be between 0 and 1")}
  
  if (zeta > 1 || zeta < 0){ 
    stop("Incorrect zeta, it must be between 0 and 1")}
  
  if (lambda > 1 || lambda < 0){
    stop("Incorrect lambda, it must be between 0 and 1")}
  
  
  parameters <- list(delta, zeta, lambda)
  names(parameters) <- c("delta", "zeta", "lambda")
  return(parameters)
}
CheckSeeds <- function(Seeds, All_genes,All_Phenotypes){
  
  Genes_Seeds_Ok <- Seeds[which(Seeds %in% All_genes)]
  Phenotype_Seeds_Ok <- Seeds[which(Seeds %in% All_Phenotypes)]
  All_seeds_ok <- c(Genes_Seeds_Ok,Phenotype_Seeds_Ok)
  All_seeds_ko <- Seeds[which(!Seeds %in% All_seeds_ok)]
  
  list_Seeds_Ok <- list(Genes_Seeds_Ok,Phenotype_Seeds_Ok)
  
  print("Seeds below do not exist in the network: ")
  print(paste(All_seeds_ko, sep=" "))
  
  if ((length(Genes_Seeds_Ok) == 0) &&  (length(Phenotype_Seeds_Ok) ==0)){
    stop("Seeds not found in our network")
  } else {
    return(list(Genes_Seeds=Genes_Seeds_Ok, Pheno_Seeds=Phenotype_Seeds_Ok))
  }
  
}

GetSeedScores <- function(geneSeeds, phenoSeeds, eta, LG, LP, tau, phi) {
  
  n <- length(geneSeeds)
  m <- length(phenoSeeds)
  
  if ((n != 0 && m != 0)) {
    
    Seed_Genes_Layer_Labeled <- paste0(rep(geneSeeds, LG), sep = "_", 
                                    rep(seq(LG), length.out = n * LG, each = n))
    Seeds_Genes_Scores <- rep(((1 - eta) * tau)/n, n)
    
    Seed_Phenos_Layer_Labeled <- paste0(rep(phenoSeeds, LP), sep = "_",
                                    rep(seq(LP),length.out = m * LP, each = m))
    Seeds_Phenos_Scores <- rep((eta * phi)/m, m)
    
  } else {
    eta <- 1
    if (n == 0) {
      Seed_Genes_Layer_Labeled <- character()
      Seeds_Genes_Scores <- numeric()
      
      Seed_Phenos_Layer_Labeled <- paste0(rep(phenoSeeds, LP), sep = "_",
                                    rep(seq(LP), length.out = m * LP, each = m))
      Seeds_Phenos_Scores <- rep((eta * phi)/m, m)
    } else {
      Seed_Genes_Layer_Labeled <- paste0(rep(geneSeeds, LG), sep = "_", 
                                    rep(seq(LG), length.out = n * LG, each = n))
      Seeds_Genes_Scores <- rep(tau/n, n)
      
      Seed_Phenos_Layer_Labeled <- character()
      Seeds_Phenos_Scores <- numeric()
    }
  }
  
  ### We prepare a data frame with the seeds.
  Seeds_Score <- data.frame(Seeds_ID = c(Seed_Genes_Layer_Labeled, 
                                         Seed_Phenos_Layer_Labeled), 
                            Score = c(Seeds_Genes_Scores, Seeds_Phenos_Scores), 
                            stringsAsFactors = FALSE)
  return(Seeds_Score)
}

CreateSupraadjacencyMatrix <- function(WholeNet, type, N, L, zeta, 
                                       is.weighted.graph) {

  Graphs <- WholeNet[which(lapply(WholeNet, `[[`, "type") == type)]
  Graphs <- lapply(Graphs, `[[`, "graph")
  
  Idem_Matrix <- Diagonal(N, x = 1)
  
  Col_Node_Names <- character()
  Row_Node_Names <- character()
  

  if (getDoParWorkers() > 1) 
    registerDoSEQ()
  cl <- makeCluster(L + 1)
  registerDoParallel(cl)
  #globalVariables("i")
  i <- 1
  SupraAdjacencyResult <- foreach(i = 1:L, .packages = c("igraph",
                                                         "Matrix")) %dopar% 
    {
      
      
      SupraAdjacencyMatrixPart <- Matrix(0, ncol = N * L, 
                                         nrow = N, sparse = TRUE)
      
      Adjacency_Layer <- as_adjacency_matrix(Graphs[[i]], attr = "weight", 
                                             sparse = TRUE)
      
      if (!is.weighted.graph) {
        Adjacency_Layer <- as_adjacency_matrix(Graphs[[i]], sparse = TRUE)
      }
      
      ## We order the matrix by the node name. This way all the matrix will have the
      ## same. Additionally we include a label with the layer number for
      ##each node name.
      Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)), 
                                         order(colnames(Adjacency_Layer))]
      Layer_Row_Col_Names <- paste(colnames(Adjacency_Layer), i, sep = "_")
      
      ## We fill the diagonal blocks with the adjacencies matrix of each layer.
      Position_ini_row <- 1
      Position_end_row <- N
      Position_ini_col <- 1 + (i - 1) * N
      Position_end_col <- N + (i - 1) * N
      SupraAdjacencyMatrixPart[(Position_ini_row:Position_end_row), 
                               (Position_ini_col:Position_end_col)] <- 
                                      (1 - zeta) * (Adjacency_Layer)
      
      # avoid division by zero for monoplex network
      L_mdfd <- L - 1
      if (L == 1) 
        L_mdfd <- 1
      
      ## We fill the off-diagonal blocks with the transition probability 
      ##among layers.
      for (j in 1:L) {
        Position_ini_col <- 1 + (j - 1) * N
        Position_end_col <- N + (j - 1) * N
        if (j != i) {
          SupraAdjacencyMatrixPart[(Position_ini_row:Position_end_row), 
                                   (Position_ini_col:Position_end_col)] <- 
                                    (zeta/(L_mdfd)) * Idem_Matrix
        }
      }
      return(list(SupraAdjacencyMatrixPart, Layer_Row_Col_Names))
    }
  
  stopCluster(cl)

  SupraAdjacencyResult <- unlist(SupraAdjacencyResult, recursive = FALSE)
  
  
  # Row-Col names are even indexed
  Col_Names <- do.call("c", SupraAdjacencyResult[seq(2, 2 * L, by = 2)])
  
  # Parallele parts of the SupraAdjacencyMatrix are odd indexed
  SupraAdjacencyMatrix <- do.call("rbind", SupraAdjacencyResult[seq(1, 2 * L,
                                                                    by = 2)])
  
  rownames(SupraAdjacencyMatrix) <- Col_Names
  colnames(SupraAdjacencyMatrix) <- Col_Names
  
  return(SupraAdjacencyMatrix)
}
CreateBipartiteMatrix <- function(WholeNet, N, M, gene_pool_nodes_sorted, 
                                  phenotype_pool_nodes_sorted, numCores) {
  Gene_Phenoivar_Network <- WholeNet[which(lapply(WholeNet, `[[`, "type") == 
                                             "bipartite")]
  Gene_Phenoivar_Network <- lapply(Gene_Phenoivar_Network, `[[`, "DF")[[1]]
  
  # Get the Subset of Gene-Phenoivar relations which have common genes in whole
  # network
  Gene_Phenoivar_Network <- Gene_Phenoivar_Network[which(
                    Gene_Phenoivar_Network$from %in% gene_pool_nodes_sorted), ]
  Gene_Phenoivar_Network <- Gene_Phenoivar_Network[which(
    Gene_Phenoivar_Network$to %in% phenotype_pool_nodes_sorted), ]
  
  Gene_Phenoivar_Network <- graph.data.frame(Gene_Phenoivar_Network, 
                                             directed = FALSE)
  
  el <- as_edgelist(Gene_Phenoivar_Network)
  value <- edge_attr(Gene_Phenoivar_Network, name = "weight")

  Bipartite_matrix <- Matrix(data = 0, nrow = N, ncol = M)
  rownames(Bipartite_matrix) <- gene_pool_nodes_sorted
  colnames(Bipartite_matrix) <- phenotype_pool_nodes_sorted
  rindx <- unlist(lapply(el[, 1], function(x) 
                      which(rownames(Bipartite_matrix) %in% x)))
  cindx <- unlist(lapply(el[, 2], function(x) 
                      which(colnames(Bipartite_matrix) %in% x)))
  
  lenRindx <- length(rindx)
  partLen <- floor(lenRindx/numCores)
  
  if (partLen == 0) {
    cat("numCores:", numCores, " - partLen: ", partLen, "\n")
    stop("Assigned numCores is greater than data length! Assign less!")
  }
  
  rindx_parts <- list()
  cindx_parts <- list()
  value_parts <- list()
  for (i in 1:numCores) {
    stInd <- (i - 1) * partLen + 1
    endInd <- i * partLen
    if (i == numCores) {
      endInd <- lenRindx
    }
    rindx_parts[[i]] <- rindx[stInd:endInd]
    cindx_parts[[i]] <- cindx[stInd:endInd]
    value_parts[[i]] <- value[stInd:endInd]
  }
  
  
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  Bipartite_matrix_result <- foreach(i = 1:numCores, 
                                     .packages = "Matrix") %dopar% 
    {
      for (j in 1:length(rindx_parts[[i]])) {
        Bipartite_matrix[rindx_parts[[i]][j], 
                         cindx_parts[[i]][j]] <- value_parts[[i]][j]
      }
      return(Bipartite_matrix)
    }
  stopCluster(cl)
  Bipartite_matrix <- Reduce("+", Bipartite_matrix_result)
  
  return(Bipartite_matrix)
}
CreateSuprabipartiteMatrix <- function(Bipartite_matrix, N, M, LG, LP) {
  SupraBipartiteMatrix <- Matrix(0, nrow = N * LG, ncol = M * LP, sparse = TRUE)
  
  Row_Node_Names <- sprintf(paste0(rep(rownames(Bipartite_matrix), LG), "_%d"), 
                            rep(seq_len(LG), each = N))
  SupraBipartiteMatrix <- do.call(rbind, replicate(LG, Bipartite_matrix,
                                                   simplify = FALSE))
  
  rownames(SupraBipartiteMatrix) <- Row_Node_Names
  
  
  Col_Node_Names <- sprintf(paste0(rep(colnames(Bipartite_matrix), LP), "_%d"), 
                            rep(seq_len(LP), each = M))
  
  SupraBipartiteMatrix <- do.call(cbind, replicate(LP, SupraBipartiteMatrix, 
                                                   simplify = FALSE))
  colnames(SupraBipartiteMatrix) <- Col_Node_Names
  return(SupraBipartiteMatrix)
}
CreateTransitionMatrix <- function(SupraBipartiteMatrix, N, M, LG, LP, lambda, 
                                     isTranspose) {
  ### TRUE = Row wise/pheno-gene, FALSE = Col Wise/gene-pheno
  if (isTranspose) {
    #### Transition Matrix for the inter-subnetworks links
    TransitionMatrix <- Matrix(0, nrow = M * LP, ncol = N * LG, sparse = TRUE)
    colnames(TransitionMatrix) <- rownames(SupraBipartiteMatrix)
    rownames(TransitionMatrix) <- colnames(SupraBipartiteMatrix)
    
    Row_Sum_Bipartite <- rowSums(SupraBipartiteMatrix, na.rm = FALSE, dims = 1, 
                                 sparseResult = FALSE)
    # row wise normalization for propability
    for (i in 1:(N * LG)) {
      if (Row_Sum_Bipartite[i] != 0) {
        TransitionMatrix[, i] <- 
                      (lambda * SupraBipartiteMatrix[i, ])/Row_Sum_Bipartite[i]
      }
    }
  } else {
    #### Transition Matrix for the inter-subnetworks links
    TransitionMatrix <- Matrix(0, nrow = N * LG, ncol = M * LP, sparse = TRUE)
    colnames(TransitionMatrix) <- colnames(SupraBipartiteMatrix)
    rownames(TransitionMatrix) <- rownames(SupraBipartiteMatrix)
    
    Col_Sum_Bipartite <- colSums(SupraBipartiteMatrix, na.rm = FALSE, dims = 1, 
                                 sparseResult = FALSE)
    
    # columnwise normalization for propability
    for (j in 1:(M * LP)) {
      if (Col_Sum_Bipartite[j] != 0) {
        TransitionMatrix[, j] <- 
                      (lambda * SupraBipartiteMatrix[, j])/Col_Sum_Bipartite[j]
      }
    }
  }
  
  
  return(TransitionMatrix)
}
CreateTransitionMatrix.gene_pheno <- function(SupraBipartiteMatrix, N, M, 
                                              LG, LP, lambda) {
  #### Transition Matrix for the inter-subnetworks links
  Transition_Gene_Phenoivar <- Matrix(0, nrow = N * LG, ncol = M * LP,
                                      sparse = TRUE)
  colnames(Transition_Gene_Phenoivar) <- colnames(SupraBipartiteMatrix)
  rownames(Transition_Gene_Phenoivar) <- rownames(SupraBipartiteMatrix)
  
  Col_Sum_Bipartite <- colSums(SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
                               sparseResult = FALSE)
  
  # columnwise normalization for propability
  for (j in 1:(M * LP)) {
    if (Col_Sum_Bipartite[j] != 0) {
      Transition_Gene_Phenoivar[, j] <- 
        (lambda * SupraBipartiteMatrix[, j])/Col_Sum_Bipartite[j]
    }
  }
  return(Transition_Gene_Phenoivar)
}
CreateTransitionMatrix.pheno_gene <- function(SupraBipartiteMatrix, N, M,
                                              LG, LP, lambda) {
  Transition_Phenoivar_Gene <- Matrix(0, nrow = M * LP, ncol = N * LG,
                                      sparse = TRUE)
  
  colnames(Transition_Phenoivar_Gene) <- rownames(SupraBipartiteMatrix)
  rownames(Transition_Phenoivar_Gene) <- colnames(SupraBipartiteMatrix)
  
  Row_Sum_Bipartite <- rowSums(SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
                               sparseResult = FALSE)
  # row wise normalization for propability
  for (i in 1:(N * LG)) {
    if (Row_Sum_Bipartite[i] != 0) {
      Transition_Phenoivar_Gene[, i] <- 
        (lambda * SupraBipartiteMatrix[i, ])/Row_Sum_Bipartite[i]
    }
  }
  return(Transition_Phenoivar_Gene)
}
CreateTransitionMultiplexNetwork <- function(SupraAdjacencyMatrix, 
                                             SupraBipartiteMatrix, Num,
                                             inputLength, lambda, numCores) {
  #### Transition Matrix for the intra-subnetworks links
  Transition_Multiplex_Network <- Matrix(0, nrow = Num * inputLength,
                                        ncol = Num * inputLength, sparse = TRUE)
  
  rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
  colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
  
  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrix, na.rm = FALSE, dims = 1, 
                               sparseResult = FALSE)
  Row_Sum_Bipartite <- rowSums(SupraBipartiteMatrix, na.rm = FALSE, dims = 1, 
                               sparseResult = FALSE)
  
  
  partLen <- floor(Num * inputLength/numCores)
  # below can only happen with toy samples
  if (partLen == 0) {
    stop("Assigned numCores is greater than data length! Assign less!")
  }
  
  parts <- vector("list", numCores)
  for (i in 1:numCores) {
    stInd <- (i - 1) * partLen + 1
    endInd <- i * partLen
    if (i == numCores) {
      endInd <- Num * inputLength
    }
    parts[[i]][["start"]] <- stInd
    parts[[i]][["end"]] <- endInd
  }
  
  cl <- makeCluster(numCores + 1)
  
  registerDoParallel(cl)
  
  Transition_Multiplex_Network_Result <- foreach(i = 1:numCores,
                                                 .packages = "Matrix") %dopar% 
    {
      for (j in parts[[i]]["start"]:parts[[i]]["end"]) {
        if (Row_Sum_Bipartite[j] != 0) {
          Transition_Multiplex_Network[, j] <- ((1 - lambda) * 
                                 SupraAdjacencyMatrix[, j])/Col_Sum_Multiplex[j]
        } else {
          Transition_Multiplex_Network[, j] <- 
            SupraAdjacencyMatrix[, j]/Col_Sum_Multiplex[j]
        }
      }
      return(Transition_Multiplex_Network)
    }
  
  Transition_Multiplex_Network <- Reduce("+", 
                                         Transition_Multiplex_Network_Result)
  stopCluster(cl)
  
  return(Transition_Multiplex_Network)
}
CreateGeneTransitionMultiplexNetwork <- function(SupraAdjacencyMatrix,
                                                 SupraBipartiteMatrix, N, LG,
                                                 lambda, numCores) {
  #### Transition Matrix for the intra-subnetworks links
  Transition_Multiplex_Network <- Matrix(0, nrow = N * LG, ncol = N * LG, 
                                         sparse = TRUE)
  
  rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
  colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
  
  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrix, na.rm = FALSE, dims = 1,
                               sparseResult = FALSE)
  Row_Sum_Bipartite <- rowSums(SupraBipartiteMatrix, na.rm = FALSE, dims = 1, 
                               sparseResult = FALSE)
  
  
  partLen <- floor(N * LG/numCores)
  # below can only happen with toy samples
  if (partLen == 0) {
    stop("Assigned numCores is greater than data length! Assign less!")
  }
  
  parts <- vector("list", numCores)
  for (i in 1:numCores) {
    stInd <- (i - 1) * partLen + 1
    endInd <- i * partLen
    if (i == numCores) {
      endInd <- N * LG
    }
    parts[[i]][["start"]] <- stInd
    parts[[i]][["end"]] <- endInd
  }
  
  cl <- makeCluster(numCores + 1)
  
  registerDoParallel(cl)
  
  Transition_Multiplex_Network_Result <- foreach(i = 1:numCores,
                                                 .packages = "Matrix") %dopar% 
    {
      for (j in parts[[i]]["start"]:parts[[i]]["end"]) {
        if (Row_Sum_Bipartite[j] != 0) {
          Transition_Multiplex_Network[, j] <- ((1 - lambda) * 
            SupraAdjacencyMatrix[, j])/Col_Sum_Multiplex[j]
        } else {
          Transition_Multiplex_Network[, j] <- 
            SupraAdjacencyMatrix[, j]/Col_Sum_Multiplex[j]
        }
      }
      return(Transition_Multiplex_Network)
    }
  
  Transition_Multiplex_Network <- Reduce("+", 
                                         Transition_Multiplex_Network_Result)
  stopCluster(cl)
  
  return(Transition_Multiplex_Network)
}
CreatePhenoTransitionMultiplexNetwork <- function(SupraAdjacencyMatrixPheno, 
                                                  SupraBipartiteMatrix, M, LP, 
                                                  lambda, numCores) {
  Transition_Multiplex_Network_Pheno <- Matrix(0, nrow = M * LP, ncol = M * LP, 
                                               sparse = TRUE)
  
  rownames(Transition_Multiplex_Network_Pheno) <-
    rownames(SupraAdjacencyMatrixPheno)
  colnames(Transition_Multiplex_Network_Pheno) <- 
    colnames(SupraAdjacencyMatrixPheno)
  
  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrixPheno, na.rm = FALSE, 
                               dims = 1, sparseResult = FALSE)
  Col_Sum_Bipartite <- colSums(SupraBipartiteMatrix, na.rm = FALSE, 
                               dims = 1, sparseResult = FALSE)
  
  partLen <- floor(M * LP/numCores)
  # below can only happen with toy samples
  if (partLen == 0) {
    stop("Assigned numCores is greater than data length! Assign less!")
  }
  
  parts <- vector("list", numCores)
  for (i in 1:numCores) {
    stInd <- (i - 1) * partLen + 1
    endInd <- i * partLen
    if (i == numCores) {
      endInd <- M * LP
    }
    parts[[i]][["start"]] <- stInd
    parts[[i]][["end"]] <- endInd
  }
  cl <- makeCluster(numCores + 1)
  
  registerDoParallel(cl)
  Transition_Multiplex_Network_Pheno_Result <- foreach(i = 1:numCores,
                                                  .packages = "Matrix") %dopar% 
    {
      for (j in parts[[i]]["start"]:parts[[i]]["end"]) {
        if (Col_Sum_Bipartite[j] != 0) {
          Transition_Multiplex_Network_Pheno[, j] <- ((1 - lambda) * 
                            SupraAdjacencyMatrixPheno[, j])/Col_Sum_Multiplex[j]
        } else {
          Transition_Multiplex_Network_Pheno[, j] <- 
            SupraAdjacencyMatrixPheno[, j]/Col_Sum_Multiplex[j]
        }
      }
      return(Transition_Multiplex_Network_Pheno)
    }
  
  Transition_Multiplex_Network_Pheno <- Reduce("+", 
                                      Transition_Multiplex_Network_Pheno_Result)
  stopCluster(cl)
  return(Transition_Multiplex_Network_Pheno)
}

rankGenes <- function(Number_Genes, Number_Layers, Results, Seeds) {
  ## We sort the score to obtain the ranking of Genes and Phenotypes.
  genes_rank <- data.frame(Gene = character(length = Number_Genes), Score = 0)
  genes_rank$Gene <- gsub("_1", "", row.names(Results)[1:Number_Genes])
  
  ## We calculate the Geometric Mean among the genes in the different layers.
  genes_rank$Score <- GeometricMean(as.vector(Results[, 1]),
                                    Number_Layers, Number_Genes)
  
  genes_rank_sort <- genes_rank[with(genes_rank, order(-Score, Gene)), ]
  
  ### We remove the seed genes from the Ranking
  genes_rank_sort_NoSeeds <- 
    genes_rank_sort[which(!genes_rank_sort$Gene %in% Seeds), ]
  
  
  genes_rank_sort_NoSeeds <- 
    genes_rank_sort_NoSeeds[, c("Gene", "Score")]
  
  return(genes_rank_sort_NoSeeds)
}
RankPhenotypes <- function(Number_Genes, Num_Gene_Layers, Number_Phenotypes, 
                           Num_Pheno_Layers, Results, Seeds) {
  
  ## rank_phenotypes
  phenotypes_rank <- data.frame(Pheno = character(length = Number_Phenotypes), 
                                Score = 0)
  phenotypes_rank$Pheno <- gsub("_1", "", 
                                row.names(Results)[(Number_Genes * 
                                Num_Gene_Layers + 1):(Number_Genes * 
                                Num_Gene_Layers + Number_Phenotypes)])
  
  phenotypes_rank$Score <- GeometricMean(as.vector(Results[, 1])[(Number_Genes * 
                                Num_Gene_Layers + 1):nrow(Results)], 
                                Num_Pheno_Layers, Number_Phenotypes)
  
  phenotypes_rank_sort <- phenotypes_rank[with(phenotypes_rank, order(-Score, 
                                                                      Pheno)), ]
  phenotypes_rank_sort_NoSeeds <- phenotypes_rank_sort[which(
                                      !phenotypes_rank_sort$Pheno %in% Seeds), ]
  
  phenotypes_rank_sort_NoSeeds$Rank <- seq(1,
                                           nrow(phenotypes_rank_sort_NoSeeds))
  phenotypes_rank_sort_NoSeeds <- phenotypes_rank_sort_NoSeeds[, c("Rank",
                                                            "Pheno", "Score")]
  
  return(phenotypes_rank_sort_NoSeeds)
}
GeometricMean <- function(Scores, L, N) {
  
  FinalScore <- numeric(length = N)
  
  for (i in seq_len(N)) {
    FinalScore[i] <- prod(Scores[seq(from = i, to = N * L, by = N)])^(1/L)
  }
  
  return(FinalScore)
}


#' @title Random Walk Restart (RWR)
#'
#' @description This method runs the random walk with restart on the provided
#'   walk matrix. It returns a data frame including ranked genes and phenotypes,
#'   and RWR scores of the genes and phenotypes. If generatePvalue is
#'   TRUE then it generates p-values along with the ranks.
#'   
#' @param walkMatrix This is the walk matrix generated by the function
#'   CreateWalkMatrix.
#' @param geneSeeds This is a vector for storing the ids of the genes that RWR starts
#' its walk. The final ranks show the proximity of the genes/phenotypes to the seed genes.
#' @param phenoSeeds This is a vector for storing the ids of the phenotypes that RWR starts
#' its walk. The final ranks show the proximity of the genes/phenotypes to the seed phenotypes. 
#' @param generatePValue If this is TRUE, it will generate the probability values for each 
#' of the gene/phenotype rankings. If it is FALSE then the function will only return the ranks of genes/phenotype.
#' @param numCores This is the number of cores used for parallel processing.
#' @param r This parameter controls the global restart probability of RWR, and it has a default value of 0.7.
#' @param eta This parameter  controls restarting of RWR either to a gene seed or phenotype seeds, 
#' higher eta means utilizing gene seeds more than phenotype seeds, and it has a default value of 0.5.
#' @param tau This is a vector that stores weights for each of the 'gene'
#'   layer in the complex gene and phenotype network. Each value of the vector
#'   corresponds to the order of the network files in the input file of CreateWalkMatrix function. The weights must
#'   sum up to the same number of gene layers. Default value gives equal weight to gene layers.
#' @param phi This is a vector that stores weights for each of the 'phenotype'
#'   layer in the complex gene and phenotype network. Each value of the vector
#'   corresponds to the order of the network files in the input file of CreateWalkMatrix function. The weights must
#'   sum up to the same number of phenotype layers. Default value gives equal weight to phenotype layers.
#'
#' @return If the parameter generatePValue is TRUE, then this function returns a
#'   data frame of ranked genes/phenotypes with p-values with three columns; Gene/Phenotype
#'   id, score, p-value. If generatePValue is FALSE, then it returns a data frame
#'   of ranked genes/phenotypes with two columns; Gene/Phenotype id, score.
#'
#' @examples
#' print("The data sets we used have not been published so the code below is unable to be run.")
#' #ranksWithPVal <- RandomWalkRestart(walkMatrix, c('gene1', 'gene2'), c(),TRUE)
#' #ranksWithPVal <- RandomWalkRestart(walkMatrix, c('gene1', 'gene2'), c(),TRUE)
#' #ranks <- RandomWalkRestart(CreateWalkMatrix('myFile.txt'),c('gene1'), 
#'  #      c('phenotype1', 'phenotype2'), FALSE)
#' #ranksWithPval <- RandomWalkRestart(CreateWalkMatrix('myFile.txt'),c('gene1'), 
#'  #      c(), TRUE, 12, 0.7, 0.6, tau =(1,0.5,1.5), phi =(1,0.5,1.5))
#' 
RandomWalkRestart <- function(walkMatrix, geneSeeds, phenoSeeds, 
                               generatePValue = TRUE, numCores = 1,
                               r = 0.7, eta = 0.5, tau = NULL, phi = NULL) {
  
  if (!is.null(tau)) {
    tau <- rep(1, walkMatrix[["LG"]])
    
    if (sum(tau)/walkMatrix[["LG"]] != 1) {
      stop("Incorrect tau, the sum of its values should be equal to the number 
      of gene layers")
    }
  }
  
  if (!is.null(phi)) {
    phi <- rep(1, walkMatrix[["LP"]])
    
    if (sum(phi)/walkMatrix[["LP"]] != 1) {
      stop("Incorrect phi, the sum of its values should be equal to the number 
      of phenotype layers")
    }
  }
  
  
  
  gene_pool_nodes_sorted <- walkMatrix[["genes"]]
  phenotype_pool_nodes_sorted <- walkMatrix[["phenotypes"]]
  
  SeedList <- CheckSeeds(c(geneSeeds, phenoSeeds), gene_pool_nodes_sorted,
                         phenotype_pool_nodes_sorted)
  Seeds_Score <- GetSeedScores(SeedList[["Genes_Seeds"]],
                               SeedList[["Pheno_Seeds"]], 
                               eta, walkMatrix[["LG"]], 
                               walkMatrix[["LP"]], tau/walkMatrix[["LG"]],
                               phi/walkMatrix[["LP"]])
  
  Threeshold <- 1e-10
  NetworkSize <- ncol(walkMatrix[["WM"]])
  
  ### We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1
  
  #We define the prox_vector(The vector we will move after the 
  #first RW iteration.
  #We start from The seed. We have to take in account that the walker with 
  #restart in some of the Seed genes, 
  #depending on the score we gave in that file).
  prox_vector <- Matrix(0, nrow = NetworkSize, ncol = 1, sparse = TRUE)
  
  prox_vector[which(colnames(walkMatrix[["WM"]]) %in% Seeds_Score[, 1])] <- 
                                                              (Seeds_Score[, 2])
  
  prox_vector <- prox_vector/sum(prox_vector)
  restart_vector <- prox_vector
  while (residue >= Threeshold) {
    old_prox_vector <- prox_vector
    prox_vector <- (1 - r) * (walkMatrix[["WM"]] %*% prox_vector) + 
                                                        r * restart_vector
    
    residue <- sqrt(sum((prox_vector - old_prox_vector)^2))
    
    iter <- iter + 1
  }
  
  RWGeneRankDF <- rankGenes(walkMatrix[["N"]], walkMatrix[["LG"]], 
                            prox_vector, SeedList[["Genes_Seeds"]])
  RWPhenoRankDF <- RankPhenotypes(walkMatrix[["N"]], walkMatrix[["LG"]], 
                                  walkMatrix[["M"]], walkMatrix[["LP"]], 
                                  prox_vector, SeedList[["Pheno_Seeds"]])
  
  if (generatePValue) {
    # Generate random seeds
    RandomSeeds <- GenerateRandomSeedVector(walkMatrix, 
                                            SeedList[["Genes_Seeds"]], 
                                            SeedList[["Pheno_Seeds"]])
    # RandomSeeds calculate random ranks walkMatrix, geneSeedsList,
    # phenoSeedsList, N, LG, LP, eta, tau, phi, r, funcs, no.cores=4
    Rand_Seed_Gene_Rank <- RandomWalkRestartBatch(walkMatrix[["WM"]], 
                                                   RandomSeeds[["gene"]], 
                                                   RandomSeeds[["phenotype"]], 
                                                   walkMatrix[["N"]], 
                                                   walkMatrix[["LG"]], 
                                                   walkMatrix[["LP"]], eta, 
                                                   tau/walkMatrix[["LG"]], 
                                                   phi/walkMatrix[["LP"]], 
                                                   r, numCores)
    
    # calculate p-values using random ranks
    dfRank <- CalculatePvalues(RWGeneRankDF, Rand_Seed_Gene_Rank, numCores)
    

    # returns cropped version will have col are gene, score, p value
    dfRankCropped <- dfRank[, c(1:2, (ncol(dfRank) - 2))]
    return(dfRankCropped)
  } else {
    rownames(RWGeneRankDF) <- NULL
    return(RWGeneRankDF)
  }
}
RandomWalkRestartSingle <- function(walkMatrix, r, Seeds_Score) {
  #We define the threshold and the number maximum of iterations for the randon
  #walker. Seeds_Score <- GetSeedScores(geneSeeds,CultSeeds, Parameters$eta, LG,
  #LC, Parameters$tau/LG, Parameters$phi/LC)
  Threeshold <- 1e-10
  NetworkSize <- ncol(walkMatrix)
  
  ### We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1
  
#We define the prox_vector(The vector we will move after the first RW iteration.
#We start from The seed. We have to take in account that the walker with restart
#in some of the Seed genes, depending on the score we gave in that file).
  prox_vector <- Matrix(0, nrow = NetworkSize, ncol = 1, sparse = TRUE)
  
  prox_vector[which(colnames(walkMatrix) %in% Seeds_Score[, 1])] <- 
                                                            (Seeds_Score[, 2])
  
  prox_vector <- prox_vector/sum(prox_vector)
  restart_vector <- prox_vector
  while (residue >= Threeshold) {
    old_prox_vector <- prox_vector
    prox_vector <- (1 - r) * (walkMatrix %*% prox_vector) + r * restart_vector
    
    residue <- sqrt(sum((prox_vector - old_prox_vector)^2))
    
    iter <- iter + 1
  }
  
  print("RWR-MH number of iteration: ")
  print(iter - 1)
  return(prox_vector)
}
GetConnectivity <- function(NetworkDF, gene_pool_nodes_sorted,
                            phenotype_pool_nodes_sorted) {
  WholeNet <- do.call("rbind", (lapply(NetworkDF, `[[`, "DF")))
  g <- graph.data.frame(WholeNet, directed = FALSE)
  A <- as_adjacency_matrix(g, sparse = TRUE, attr = "weight")
  Degree = apply(A, 2, sum)
  Connectivity <- data.frame(Node = as.character(A@Dimnames[[1]]), 
                             Degree = Degree, row.names = NULL)
  Connectivity <- Connectivity[order(Connectivity$Degree, decreasing = TRUE), ]
  Connectivity$Node <- as.character(Connectivity$Node)
  GeneConnectivity <- Connectivity[which(Connectivity$Node %in% 
                                           gene_pool_nodes_sorted), ]
  PhenoConnectivity <- Connectivity[which(Connectivity$Node %in% 
                                            phenotype_pool_nodes_sorted), ]
  return(list(gene = GeneConnectivity, pheno = PhenoConnectivity))
}


#' @title Create Walk Matrix
#'
#' @description Generates a Walk Matrix (Transition Matrix) from Gene and Phenotype networks for RWR.
#'
#' @param inputFileName The name of the text file that contains the names of
# gene and phenotype network files. This text file is made
# from two columns that are tab-separated. The first row is header with
# 'type' and 'file_name'. Network files can be either gene, phenotype, or bipartite files.
# Each network file is composed of three columns 'from', 'to', and 'weight' columns
# which are all tab-separated. For gene/phenotype network layers, the 'from' and
# 'to' columns will have genes and phenotypes. The weight
# column will have the weight of the interactions of the genes/phenotypes,
# for unweighted network all weights must have a value of 1. For a bipartite
# network file, the 'from' column must have genes, the 'to' column must have phenotypes, the
# meaning and usage of weight is similar to the gene and phenotype layers. For
# more information on the file formatting, please check the vignette. 
#' @param numCores This is the number of cores used for parallel processing.
#'
#' @param delta This is the probability of jumping between gene layers. High delta means 
#' RWR  is high likely to jump to other layers in gene multiplex netwrok. It has a default value of 0.5.
#' @param zeta This is the probability of jumping between phenotype layers. High zeta means 
#' RWR  is high likely to jump to other layers in phenotype multiplex netwrok. It has a default value of 0.5.
#' @param lambda This is the probability of jumpting between gene and phenotype multiplex networks. High lambda
#' means RWR is more likely to exploit the bipartite relation. It has a default value of 0.5.
#'
                           
#' @return This returns a list containing the walk matrix, a sorted list of gene ids, a sorted list of phenotype ids,
#'  the connectivity degree of the genes, the connectivity degree of the phenotypes, the number of gene layers, 
#'  the number of phenotype layers, the number of genes and the number of phenotypes in the final complex network.
#'   
#' @export
#'
#' @examples
#' print("The data sets we used have not been published so the code below is unable to be run")
#' #CreateWalkMatrix('myInput.txt')
#' #CreateWalkMatrix('file.txt', detectCores(), 0.4, 0.7, 0.9)
#' 
#' 
CreateWalkMatrix <- function(inputFileName, numCores = 1, delta = 0.5, 
                             zeta = 0.5, lambda = 0.5) {
  
  Settings <- ReadSettings(inputFileName)
  Parameters <- AddParameters(delta, zeta, lambda)
  FilesDF <- Settings$FilesDF
  LG <- Settings$LG
  LP <- Settings$LP
  
  FullNet <- ReadNetworkLayers(FilesDF)

  gene_pool_nodes_sorted <- FullNet$gene_pool_nodes_sorted
  phenotype_pool_nodes_sorted <- FullNet$phenotype_pool_nodes_sorted
  FullNet <- FullNet$NetworkLayers
  
  N <- length(gene_pool_nodes_sorted)
  M <- length(phenotype_pool_nodes_sorted)

  SupraAdjacencyMatrix <- CreateSupraadjacencyMatrix(FullNet, "gene", N, LG, 
                                                       Parameters$zeta, TRUE)

  SupraAdjacencyMatrixPheno <- CreateSupraadjacencyMatrix(FullNet, "phenotype", 
                                                          M, LP, 
                                                          Parameters$delta,
                                                          TRUE)
  BipartiteMatrix <- CreateBipartiteMatrix(FullNet, N, M, 
                                          gene_pool_nodes_sorted, 
                                          phenotype_pool_nodes_sorted, numCores)
  
  ## We expand the biparite graph to fit the multiplex dimensions.  The biparti
  ## matrix has now NL x MK The genes in all the layers have to point to the
  ## phenotypes in all layers
  
  SupraBipartiteMatrix <- CreateSuprabipartiteMatrix(BipartiteMatrix, N, M, LG, 
                                                       LP)
  # Parameters$lambda)
  Transition_Gene_Phenoivar <- CreateTransitionMatrix(SupraBipartiteMatrix, N, 
                                                      M, LG, LP, 
                                                      Parameters$lambda, FALSE)
  Transition_Phenoivar_Gene <- CreateTransitionMatrix(SupraBipartiteMatrix, N, 
                                                      M, LG, LP, 
                                                      Parameters$lambda, TRUE)
  Gene_Transition_Multiplex_Network <- CreateGeneTransitionMultiplexNetwork(
                                                          SupraAdjacencyMatrix, 
                                                          SupraBipartiteMatrix, 
                                                          N, LG, 
                                                          Parameters$lambda, 
                                                          numCores)

  # CREATE TRANSITION MULTIPLEX NETWORK FOR CULTIVARS t1 <- Sys.time()
  Pheno_Transition_Multiplex_Network <- CreatePhenoTransitionMultiplexNetwork(
                                                      SupraAdjacencyMatrixPheno, 
                                                      SupraBipartiteMatrix, M, 
                                                      LP, Parameters$lambda, 
                                                      numCores)

    Multiplex_Heterogeneous_Matrix <- rbind(cbind(
                                            Gene_Transition_Multiplex_Network, 
                                            Transition_Gene_Phenoivar), 
                                            cbind(Transition_Phenoivar_Gene, 
                                            Pheno_Transition_Multiplex_Network))
  
  # Extract candidate genes for further Random Walks on this WM
  GenePhenoDF <- lapply(FullNet[which(lapply(FullNet, `[[`, "type") == 
                                            "bipartite")], `[[`, "DF")[[1]]
  CandidateGenes <- unique(GenePhenoDF$from)
  Connectivity <- GetConnectivity(FullNet, gene_pool_nodes_sorted, 
                                  phenotype_pool_nodes_sorted)
  WM <- list(WM = Multiplex_Heterogeneous_Matrix, 
             genes = gene_pool_nodes_sorted, 
             phenotypes = phenotype_pool_nodes_sorted, 
             gene_connectivity = Connectivity[["gene"]], 
             phenotype_connectivity = Connectivity[["pheno"]], 
             LG = LG, LP = LP, N = N, M = M)
  
  registerDoSEQ()
  return(WM)
}



AssignGroupToConnectivityDF <- function(ConnectivityDF, no.groups) {
  chunk.size <- ceiling(nrow(ConnectivityDF)/no.groups)
  groups <- rep(1:no.groups, each = chunk.size, length.out =
                  nrow(ConnectivityDF))
  ConnectivityDF$Group <- groups
  ConnectivityDF
}
GenerateRandomSeeds <- function(Seeds, ConnectivityDF, S = 1000, no.groups = 10, 
                                  replace_bool = FALSE) {
  seed.set.size <- length(Seeds)
  sample_size <- ceiling((S/no.groups) * seed.set.size)
  #set.seed(1)

  # Stratified Sample 'sample_size' nodes from each group as Random Seeds
  ConnectivityDF <- ConnectivityDF[which(!ConnectivityDF$Node %in% Seeds), ] 
  ConnectivityDF <- dplyr::group_by(ConnectivityDF, ConnectivityDF$Group)
  RandomSeeds <-  dplyr::sample_n(ConnectivityDF, sample_size, 
                                  replace = replace_bool)
  
  #We are creating 'seed.set.size' length seed sets by taking 'nodes' from each
  #group To this end, we determine a split order and sort the RandomSeeds DF wrt
  #to this 'order' column
  order_vec <- vector()
  for (i in 1:no.groups) {
    order_vec <- c(order_vec, seq(from = i, to = S * seed.set.size,
                                  by = no.groups))
  }
  RandomSeeds$Order <- as.factor(order_vec)
  RandomSeeds <- RandomSeeds[order(RandomSeeds$Order), ]
  
  # split the sorted DF into 'seed.set.size' length vectors
  RandomSeeds <- split(RandomSeeds$Node, 
                       (seq(nrow(RandomSeeds)) - 1)%/%seed.set.size)
  
  RandomSeeds
}

GenerateRandomSeedVector <- function(WM, geneSeeds, phenoSeeds, S = 1000, 
                                     no.groups.gene = 10, 
                                     no.groups.pheno = 5) {
  
  if (length(geneSeeds) == 0 && length(phenoSeeds) == 0) 
    stop("No seeds provided!")
  GeneConnectivity <- AssignGroupToConnectivityDF(WM[["gene_connectivity"]], 
                                                     no.groups = no.groups.gene)
  PhenoConnectivity<-AssignGroupToConnectivityDF(WM[["phenotype_connectivity"]], 
                                                    no.groups = no.groups.pheno)

  RandomgeneSeeds <- list()
  if (length(geneSeeds) != 0) {
    RandomgeneSeeds <- GenerateRandomSeeds(geneSeeds, GeneConnectivity, S, 
                                             no.groups.gene, TRUE)
    if (length(geneSeeds) != 1 && any(duplicated(RandomgeneSeeds[1:S]))) 
      warning("WARN: Duplicated random 'gene' seeds generated!")
  }
  
  RandomphenoSeeds <- list()
  if (length(phenoSeeds) != 0) {
    RandomphenoSeeds <- GenerateRandomSeeds(phenoSeeds, PhenoConnectivity, 
                                              S, no.groups.pheno, TRUE)
    if (length(phenoSeeds) != 1 && any(duplicated(RandomphenoSeeds[1:S]))) 
      warning("WARN: Duplicated random 'phenotype' seeds generated!")
  }
  
  return(list(gene = RandomgeneSeeds, phenotype = RandomphenoSeeds))
  
}


CalculatePvalues <- function(RWGeneRanks, Rand_Seed_Gene_Rank, no.cores) {
  # t <- Sys.time()
  S <- ncol(Rand_Seed_Gene_Rank)/3
  cl <- makeCluster(no.cores)
  registerDoParallel(cl)
  dfRanks <- foreach(i = 1:(S)) %dopar% {
    # traverse all gene names and get their ranks in random run result
    rand_ranks <- sapply(RWGeneRanks$Gene, function(gene) {
      which(Rand_Seed_Gene_Rank[, 2 + 3 * (i - 1)] %in% gene)
    })
    
    #if genes are in the random seed set for this run then there are no rank for
    #them
    idx <- !(sapply(rand_ranks, length))
    rand_ranks[idx] <- NA
    dfRank <- unname(unlist(rand_ranks))
    dfRank
  }
  stopCluster(cl)
  # getDoParWorkers() create dfRanks DF from list output of random seed ranks
  dfRanks <- suppressMessages(as.data.frame(bind_cols(dfRanks)))
  
  # add the Gene names as the first column
  dfRanks <- cbind(RWGeneRanks$Gene, RWGeneRanks$Score, dfRanks, 
                   stringsAsFactors = FALSE)
  colnames(dfRanks)[1:2] <- c("Gene", "Score")
  
  rank.offset <- seq(10, 100, by = 10)
  # create offset rank columns by adding offset vector for comparison
  dfRanks <- cbind(dfRanks, replicate(length(rank.offset), 
                      as.numeric(rownames(dfRanks))) + 
                      t(replicate(nrow(dfRanks), rank.offset)))
  colnames(dfRanks)[(S + 3):(S + length(rank.offset) + 2)] <- paste0("Rank", 
                                                                    rank.offset)
  
  # calculate P values by comparing (random seed rank + offset value) vs (actual
  # seed rank)
  for (i in 1:length(rank.offset)) {
    dfRanks["P_value"] <- base::rowMeans(dfRanks[, 3:(S+2)] < 
                                                (dfRanks[,S+2+i]), na.rm = TRUE) 
  }
  
  
  # Calculate Median and Average ranks of genes for Random Seeds
  dfRanks$Med <- apply(dfRanks[, 3:(2 + S)], 1, median)
  dfRanks$Ave <- rowMeans(dfRanks[, 3:(2 + S)])
  dfRanks
}

RandomWalkRestartBatch <- function(walkMatrix, geneSeedsList, phenoSeedsList, 
                                    N, LG, LP, eta, tau, phi, r, no.cores = 4) {
  cl <- makeCluster(no.cores)
  registerDoParallel(cl)
  seedsLength <- ifelse(length(geneSeedsList) != 0, length(geneSeedsList), 
                        length(phenoSeedsList))
  funcs <- c("GetSeedScores", "RandomWalkRestartSingle", "rankGenes", 
             "GeometricMean")
  #globalVariables("i")
  i <- 1
  Rand_Seed_Gene_Rank <- foreach(i = 1:seedsLength, .combine = cbind, .export =
                                   funcs, .packages = c("Matrix")) %dopar% {
      

      if (length(geneSeedsList) != 0 && length(phenoSeedsList) != 0) {
          Seeds_Score <- GetSeedScores(geneSeedsList[[i]], phenoSeedsList[[i]], 
                                                        eta, LG, LP, tau, phi)
      } else if (length(geneSeedsList) != 0) {
          Seeds_Score <- GetSeedScores(geneSeedsList[[i]], vector(), eta, LG, 
                                                                LP, tau, phi)
      } else {
          Seeds_Score <- GetSeedScores(vector(), phenoSeedsList[[i]], eta, LG, 
                                                                LP, tau, phi)
      }
                                   
      Rand_Seed_Res <- RandomWalkRestartSingle(walkMatrix, r, Seeds_Score)
      Rand_Seed_Gene_Rank <- rankGenes(N, LG, Rand_Seed_Res, 
                                       ifelse(length(geneSeedsList) !=  
                                            0, geneSeedsList[[i]], vector()))
      return(Rand_Seed_Gene_Rank)
    }
  stopCluster(cl)

  return(Rand_Seed_Gene_Rank)
}
