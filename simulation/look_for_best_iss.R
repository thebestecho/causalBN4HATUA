

library(bnlearn)
library(MCMCpack)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)

### Load Functions
RBIND <- function(datalist, keep.rownames = TRUE) {
  Len <- sapply(datalist, ncol)
  if (all(diff(Len) == 0)) {
    temp <- names(datalist[[1]])
    if (all(sapply(datalist, function(x) names(x) %in% temp))) tryme <- "basic"
    else tryme <- "complex"
  } 
  else tryme <- "complex"
  almost <- switch(
    tryme,
    basic = { do.call("rbind", datalist) },
    complex = {
      Names <- unique(unlist(lapply(datalist, names)))
      NROWS <- c(0, cumsum(sapply(datalist, nrow)))
      NROWS <- paste(NROWS[-length(NROWS)]+1, NROWS[-1], sep=":")
      out <- lapply(1:length(datalist), function(x) {
        emptyMat <- matrix(NA, nrow = nrow(datalist[[x]]), ncol = length(Names))
        colnames(emptyMat) <- Names
        emptyMat[, match(names(datalist[[x]]), 
                         colnames(emptyMat))] <- as.matrix(datalist[[x]])
        emptyMat
      })
      do.call("rbind", out)
    })
  Final <- as.data.frame(almost, row.names = 1:nrow(almost))
  Final <- data.frame(lapply(Final, function(x) type.convert(as.character(x))))
  if (isTRUE(keep.rownames)) {
    row.names(Final) <- make.unique(unlist(lapply(datalist, row.names)))
  } 
  Final
}

## generate a random network structure with specified number of nodes and maximum number of parents per node
randomBN_structure <- function(nodes, maxin = 3, ...) {
  
  # ic-dag: Ide's and Cozman's Generating Multi-connected DAGs algorithm
  if(is.numeric(nodes)){
    numnodes <- nodes
    nodeNames <- as.character(1:numnodes)
    ranBN <- random.graph(nodeNames, method = "ic-dag", 
                          max.in.degree = maxin, ... = ...)
  } else {
    numnodes <- length(nodes)
    nodeNames <- nodes
    ranBN <- random.graph(nodes, method = "ic-dag", max.in.degree = maxin, ... = ...)
  }
  
  return(ranBN)
  
}


## assign parameters to the network with varying levels for discrete nodes
randomBN_paras <- function(ranBN, levelsVec = levelsVec, alphaparent = 0.5, alphanone = 5) {
  
  nodeNames <- nodes(ranBN)
  
  # get the CPTs for each node
  CPTs <- vector("list", length = length(nodeNames))
  names(CPTs) <- nodeNames
  
  for(node in nodeNames){
    indegree <- in.degree(ranBN, node) # Ensure in.degree is properly defined or replaced
    nlevels <- levelsVec[which(nodeNames == node)]
    
    if(indegree == 0) { 
      
      # No parents
      alpha_vector <- rep(alphanone, nlevels)
      CPTs[[node]] <- rdirichlet(1,alpha_vector)
      
    } else { 
      
      # Has parents
      parentNodes <- parents(ranBN, node) 
      
      # Ensure parents is properly defined or replaced
      if (length(parentNodes) == 0) {
        
        parentStatesComb <- 1
        
      } else {
        parentLevelsIndices <- match(parentNodes, nodeNames)
        parentStatesComb <- prod(levelsVec[parentLevelsIndices])
      }
      
      if(parentStatesComb <= 0) {
        stop("Invalid parent states combination calculated.")
      } else {
        
        # Initializing an array for the CPT with correct dimensions
        alpha_vector <- rep(alphaparent, nlevels)
        tempdirichlet <- rdirichlet(parentStatesComb,alpha_vector) ### this is the key
        transdirichlet <- t(tempdirichlet)
        dimsdirichlet <- c(nlevels, levelsVec[parentLevelsIndices])
        dim(transdirichlet) <- dimsdirichlet
        CPTs[[node]] <- transdirichlet
      }
    }
  }
  
  # and now make the parameterised BN
  ranBNparam <- custom.fit(ranBN, dist = CPTs)
  return(ranBNparam)
}



## find the iss that gives the highest score, range 1-100 with interval 1
best_iss <- function(seed_structure, seed_data, levelsVec, n_vars, npoint, maxin = 3){
  
  ### Step 1 - generate a random network structure
  set.seed(seed_structure)
  ranBN <- randomBN_structure (nodes = n_vars, maxin = maxin)
  
  ### Step 2 - parameterize the network with specified levels
  BNinfo <- randomBN_paras(ranBN = ranBN, levelsVec = levelsVec, alphaparent = 0.5, alphanone = 5)
  
  ### Step 3 - simulate the data set based on the parameterized network
  set.seed(seed_data)
  simulated_data <- rbn(BNinfo, npoint) %>% sapply(as.integer) %>% as.data.frame()
  fac_cols <- colnames(simulated_data) # Covert to factor
  simulated_data <- simulated_data %>% mutate(across(all_of(fac_cols), as.factor))
  
  ### Step 4 - score the networks structure with simulated data set using varying iss values
  ### from 1 to 100 with interval 1
  results <- lapply(seq(1,100,length.out=100),function(x){
    iss <- x
    scores <- score(ranBN, simulated_data, type = "bde", iss=x)
    cbind(scores,iss) %>% as.data.frame()
  })
  
  results <- results %>% RBIND(keep.rownames = TRUE)
  
  ### extract the iss that gives the highest score
  best <- results$iss[results$scores == max(results$scores)]
  
  return(best)
  
}


## 100 different data sets for different original network structure
## each random network structure has one single data set
myfun <- function(structure_seed, data_seed, levelsVec, n_vars, npoint, maxin, data_type){
s <- data_seed
var1 <- parallel::parLapply(cluster,s:(s+99),function(x){

  best_iss(seed_structure=x, seed_data=x, levelsVec=levelsVec, n_vars=n_vars, npoint=npoint, maxin = maxin)
  })
  
var1 <- unlist(var1) %>% as.data.frame()
write.csv(var1,paste(data_type,".csv",sep = ""))
}

## update when used for different subsets
levelsVec <- c(9, 2, 4, 4, 2,
               6, 2, 2, 4, 2,
               2, 2, 3, 2, 3,
               2, 2, 2, 4, 2,
               2, 2, 2, 2, 2,
               2, 2, 2)

library(parallel)
library(doSNOW)
cluster <- parallel::makeCluster(48)

parallel::clusterExport(cluster,"RBIND")
parallel::clusterExport(cluster,"randomBN_structure")
parallel::clusterExport(cluster,"randomBN_paras")
parallel::clusterExport(cluster,"best_iss")
parallel::clusterExport(cluster,"myfun")
parallel::clusterExport(cluster,"levelsVec")
parallel::clusterEvalQ(cluster,library(bnlearn))
parallel::clusterEvalQ(cluster,library(MCMCpack))
parallel::clusterEvalQ(cluster,library(dplyr))
parallel::clusterEvalQ(cluster,library(tidyr))
parallel::clusterEvalQ(cluster,library(purrr))
parallel::clusterEvalQ(cluster,library(magrittr))

registerDoSNOW(cluster)

args <- commandArgs(trailingOnly = TRUE)


### overall with site               
myfun(structure_seed=10001, data_seed=10001, levelsVec=levelsVec, n_vars=28, npoint=2007, maxin=3, data_type="overall_with")

levelsVec <- c(9, 2, 4, 4, 2,
               6, 2, 2, 4, 2,
               2, 2, 3, 2, 3,
               2, 2, 2, 4, 2,
               2, 2, 2, 2, 2,
               2, 2)

parallel::clusterExport(cluster,"levelsVec")

### female with site               
myfun(structure_seed=10001, data_seed=10001, levelsVec=levelsVec, n_vars=27, npoint=1722, maxin=3, data_type="female_with")

### male with site               
myfun(structure_seed=10001, data_seed=10001, levelsVec=levelsVec, n_vars=27, npoint=285, maxin=3, data_type="male_with")

levelsVec <- c(2, 4, 4, 2,
               6, 2, 2, 4, 2,
               2, 2, 3, 2, 3,
               2, 2, 2, 4, 2,
               2, 2, 2, 2, 2,
               2, 2, 2)

parallel::clusterExport(cluster,"levelsVec")

### overall with site               
myfun(structure_seed=10001, data_seed=10001, levelsVec=levelsVec, n_vars=27, npoint=2007, maxin=3, data_type="overall_without")

levelsVec <- c(2, 4, 4, 2,
               6, 2, 2, 4, 2,
               2, 2, 3, 2, 3,
               2, 2, 2, 4, 2,
               2, 2, 2, 2, 2,
               2, 2)

parallel::clusterExport(cluster,"levelsVec")

### female with site               
myfun(structure_seed=10001, data_seed=10001, levelsVec=levelsVec, n_vars=26, npoint=1722, maxin=3, data_type="female_without")

### male with site               
myfun(structure_seed=10001, data_seed=10001, levelsVec=levelsVec, n_vars=26, npoint=285, maxin=3, data_type="male_without")
              

# stop cluster and remove clients
parallel::stopCluster(cluster); print("Cluster stopped.")

# insert serial backend, otherwise error in repetitive tasks
registerDoSEQ()
remove(cluster)

