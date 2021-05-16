nullcomm_FD <- function(com,type="all"){
  if(type == "col-wise"){  # Relocate abundances within a column
    for (i in 1:ncol(com)) {
      com[,i] <- sample(com[,i], size = length(com[,i]), replace = F)
    }
  }
  if(type == "row-wise"){ # Relocate abundances within a row
    for (i in 1:nrow(com)) {
      com[i,] <- sample(com[i,], size = length(com[i,]), replace = F)
    }
  }
  if(type == "all"){ # Relocate abundances randomly within the matrix
    for (i in 1:ncol(com)) {
      com[,i] <- sample(com[,i], size = length(com[,i]), replace = F)
    }
    for (i in 1:nrow(com)) {
      com[i,] <- sample(com[i,], size = length(com[i,]), replace = F)
    }
    
  }
  return(com)
  
}
