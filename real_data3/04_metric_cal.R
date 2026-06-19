rm(list = ls())
library(igraph)
library(pROC)

source("functions/data_process.R")
source("real_data3/01_functions.R")


friendship_env <- new.env()
substances_env <- new.env()
lifestyle_env <- new.env()
selections_env <- new.env()
# 
load("datas/Glasgow-friendship.RData", envir = friendship_env)
load("datas/Glasgow-substances.RData", envir = substances_env)
load("datas/Glasgow-lifestyle.RData", envir = lifestyle_env)
load("datas/Glasgow-selections.RData", envir = selections_env)
selections = selections_env$selection129
selections <- selections_env$selection129

n <- sum(selections)

F1 <- friendship_env$friendship.1[
  selections,
  selections
]

A <- matrix(0, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    if (F1[i,j] + F1[j,i] > 0)
      A[i,j] <- 1
  }
}

Y1 <- matrix(
  substances_env$alcohol[selections,1],
  n,1
)

MissY1 <- is.na(Y1) * 1
Y1[is.na(Y1)] <- 1

Y2 <- matrix(
  substances_env$tobacco[selections,1],
  n,1
)

Y3 <- cbind(
  substances_env$cannabis[selections,1],
  lifestyle_env$leisure1[selections,]
)

MissY3 <- is.na(Y3) * 1
Y3[is.na(Y3)] <- 1

Y4 <- lifestyle_env$music1[selections,]


g <- graph_from_adjacency_matrix(
  A,
  mode = "undirected"
)

global_trans <- transitivity(
  g,
  type = "global"
)

cat("Observed transitivity =", global_trans, "\n")


# load result


trans3 <- numeric(100)
trans5 <- numeric(100)
trans7 <- numeric(100)

for (i in 1:100) {
  
  trans3_post <- numeric(1000)
  trans5_post <- numeric(1000)
  trans7_post <- numeric(1000)
  
  for (j in 1:1000) {
    
    ## k = 5
    
    alphai <- res[[i]]$alpha[j,]
    Zi <- res[[i]]$Z[[j]]
    
    Ai <- get_A(alphai,Zi)
    
    gi <- graph_from_adjacency_matrix(
      Ai,
      mode = "undirected"
    )
    
    trans5_post[j] <- transitivity(
      gi,
      type = "global"
    )
    
    ## k = 3
    
    alphai <- res3[[i]]$alpha[j,]
    Zi <- res3[[i]]$Z[j,,]
    
    Ai <- get_A(alphai,Zi)
    
    gi <- graph_from_adjacency_matrix(
      Ai,
      mode = "undirected"
    )
    
    trans3_post[j] <- transitivity(
      gi,
      type = "global"
    )
    
    ## k = 7
    
    alphai <- res7[[i]]$alpha[j,]
    Zi <- res7[[i]]$Z[j,,]
    
    Ai <- get_A(alphai,Zi)
    
    gi <- graph_from_adjacency_matrix(
      Ai,
      mode = "undirected"
    )
    
    trans7_post[j] <- transitivity(
      gi,
      type = "global"
    )
    
  }
  
  trans3[i] <- mean(trans3_post)
  trans5[i] <- mean(trans5_post)
  trans7[i] <- mean(trans7_post)
}


trans_error <- cbind(
  k3 = trans3,
  k5 = trans5,
  k7 = trans7
) - global_trans

cat("\nMSE of transitivity:\n")
print(colMeans(trans_error^2))

cat("\nSD:\n")
print(apply(trans_error,2,sd))


Loc <- upper.tri(A)

AUC3 <- numeric(100)
AUC5 <- numeric(100)
AUC7 <- numeric(100)

for (i in 1:100) {
  
  AUC3_post <- array(0,c(1000,n,n))
  AUC5_post <- array(0,c(1000,n,n))
  AUC7_post <- array(0,c(1000,n,n))
  
  for (j in 1:1000) {
    
    alphai <- res[[i]]$alpha[j,]
    Zi <- res[[i]]$Z[[j]]
    
    AUC5_post[j,,] <- plogis(
      Zi %*% t(Zi) +
        matrix(1,n,1) %*%
        matrix(alphai,1,n) +
        matrix(alphai,n,1) %*%
        matrix(1,1,n)
    )
    
    alphai <- res3[[i]]$alpha[j,]
    Zi <- res3[[i]]$Z[j,,]
    
    AUC3_post[j,,] <- plogis(
      Zi %*% t(Zi) +
        matrix(1,n,1) %*%
        matrix(alphai,1,n) +
        matrix(alphai,n,1) %*%
        matrix(1,1,n)
    )
    
    alphai <- res7[[i]]$alpha[j,]
    Zi <- res7[[i]]$Z[j,,]
    
    AUC7_post[j,,] <- plogis(
      Zi %*% t(Zi) +
        matrix(1,n,1) %*%
        matrix(alphai,1,n) +
        matrix(alphai,n,1) %*%
        matrix(1,1,n)
    )
    
  }
  
  AUC5[i] <- auc(
    roc(
      A[Loc],
      apply(AUC5_post,c(2,3),mean)[Loc],
      quiet = TRUE
    )
  )
  
  AUC3[i] <- auc(
    roc(
      A[Loc],
      apply(AUC3_post,c(2,3),mean)[Loc],
      quiet = TRUE
    )
  )
  
  AUC7[i] <- auc(
    roc(
      A[Loc],
      apply(AUC7_post,c(2,3),mean)[Loc],
      quiet = TRUE
    )
  )
}

AUC <- cbind(
  k3 = AUC3,
  k5 = AUC5,
  k7 = AUC7
)

cat("\nMean AUC:\n")
print(colMeans(AUC))

cat("\nSD AUC:\n")
print(apply(AUC,2,sd))



# load result
trans_coss <- numeric(100)
HH <- numeric(100)

for (i in 1:100) {
  
  trans_post <- numeric(1000)
  
  for (j in 1:1000) {
    
    alphai <- res[[i]]$alpha[j,]
    Zi <- res[[i]]$Z[j,,]
    
    Ai <- get_A(alphai,Zi)
    
    gi <- graph_from_adjacency_matrix(
      Ai,
      mode = "undirected"
    )
    
    trans_post[j] <- transitivity(
      gi,
      type = "global"
    )
    
  }
  
  
  trans_coss[i] <- mean(trans_post)
}

cat("\nCOSS MSE:\n")
print(mean((trans_coss - global_trans)^2))

cat("\nCOSS SD:\n")
print(sd((trans_coss - global_trans)^2))
