source("real_data3/01_functions.R")
nrun = 15000
burn = 10000
thin = 5
source("functions/data_process.R")
H_all = matrix(0, replation, nrun)
for(i in 1:replation){
  H_all[i, ] = res[[i]]$H
}
thining = seq(burn + thin, nrun, thin)
H_thin = H_all[, thining]
H_hat = apply(H_thin, 1, getmode)
table(H_hat)

# load result before
H_all = matrix(0, replation, nrun)
for(i in 1:replation){
  H_all[i, ] = res[[i]]$H
}
thining = seq(burn + thin, nrun, thin)
H_thin = H_all[, thining]
H_hat = apply(H_thin, 1, getmode)


get_model_selection_data3 <- function(res, A, Y1, Y2, Y3, Y4, MissY1, MissY3,
                                      max_dim = 10, n_rep = 100) {
  
  AIC  <- matrix(NA, max_dim, n_rep)
  BIC  <- matrix(NA, max_dim, n_rep)
  DIC  <- matrix(NA, max_dim, n_rep)
  WAIC <- matrix(NA, max_dim, n_rep)
  
  for (i in 1:n_rep) {
    
    cat("Replication:", i, "\n")
    
    for (m in 1:max_dim) {
      
      fit <- res[[i]][[m]]
      
      # posterior means for AIC/BIC
      alpha  <- apply(fit$alpha, 2, mean)
      b11    <- apply(fit$b11, 2, mean)
      b12    <- apply(fit$b12, 2, mean)
      b13    <- apply(fit$b13, 2, mean)
      b14    <- apply(fit$b14, 2, mean)
      b21    <- apply(fit$b21, 2, mean)
      b22    <- apply(fit$b22, 2, mean)
      b31    <- apply(fit$b31, 2, mean)
      b32    <- apply(fit$b32, 2, mean)
      b33    <- apply(fit$b33, 2, mean)
      gamma4 <- apply(fit$gamma4, 2, mean)
      
      Z  <- apply(fit$Z,  c(2,3), mean)
      B1 <- apply(fit$B1, c(2,3), mean)
      B2 <- apply(fit$B2, c(2,3), mean)
      B3 <- apply(fit$B3, c(2,3), mean)
      B4 <- apply(fit$B4, c(2,3), mean)
      
      # AIC
      AIC[m,i] <- get_AIC_data3(
        A,Y1,Y2,Y3,Y4,MissY1,MissY3,
        b11,b12,b13,b14,b21,b22,
        b31,b32,b33,gamma4,
        B1,B2,B3,B4,alpha,Z
      )
      
      # BIC
      BIC[m,i] <- get_BIC_data3(
        A,Y1,Y2,Y3,Y4,MissY1,MissY3,
        b11,b12,b13,b14,b21,b22,
        b31,b32,b33,gamma4,
        B1,B2,B3,B4,alpha,Z
      )
      
      # DIC
      DIC[m,i] <- get_DIC_data3(
        A,Y1,Y2,Y3,Y4,MissY1,MissY3,
        fit$b11, fit$b12, fit$b13, fit$b14,
        fit$b21, fit$b22,
        fit$b31, fit$b32, fit$b33,
        fit$gamma4,
        fit$B1, fit$B2, fit$B3, fit$B4,
        fit$alpha, fit$Z
      )
      
      # WAIC
      WAIC[m,i] <- get_WAIC_data3(
        A,Y1,Y2,Y3,Y4,MissY1,MissY3,
        fit$b11, fit$b12, fit$b13, fit$b14,
        fit$b21, fit$b22,
        fit$b31, fit$b32, fit$b33,
        fit$gamma4,
        fit$B1, fit$B2, fit$B3, fit$B4,
        fit$alpha, fit$Z
      )
    }
  }
  
  # selected dimensions
  AIC_min  <- apply(AIC,  2, which.min)
  BIC_min  <- apply(BIC,  2, which.min)
  DIC_min  <- apply(DIC,  2, which.min)
  WAIC_min <- apply(WAIC, 2, which.min)
  
  selection_table <- list(
    AIC  = table(AIC_min),
    BIC  = table(BIC_min),
    DIC  = table(DIC_min),
    WAIC = table(WAIC_min)
  )
  
  return(list(
    AIC = AIC,
    BIC = BIC,
    DIC = DIC,
    WAIC = WAIC,
    AIC_min = AIC_min,
    BIC_min = BIC_min,
    DIC_min = DIC_min,
    WAIC_min = WAIC_min,
    selection_table = selection_table
  ))
}

# load result before
model_selection <- get_model_selection_data3(
  res,
  A, Y1, Y2, Y3, Y4,
  MissY1, MissY3
)

model_selection$selection_table$AIC
model_selection$selection_table$BIC
model_selection$selection_table$DIC
model_selection$selection_table$WAIC

# load result before

replation = 100
loglikehood_total = matrix(0,10,replation)
for (i in 1:replation) {
  loglikehood_total[,i] =  apply(res[[i]], 1,mean)
}
loglikehood_total_max = apply(loglikehood_total, 2, which.max)
print(table(loglikehood_total_max))
loglikehood_total1 = loglikehood_total[1:7,]
loglikehood_best = apply(loglikehood_total1, 2, max)
loglikehood_sd = apply(loglikehood_total1, 2, sd)
for(i in 1:replation){
  loglikehood_total1[(loglikehood_total1[,i] < loglikehood_best[i] - loglikehood_sd[i]),i] = 0
}
find_min_nonzero_indices <- function(mat) {
  apply(mat, 2, function(col) {
    nonzero_indices <- which(col != 0)
    if (length(nonzero_indices) > 0) {
      min(nonzero_indices)
    } else {
      NA  # SS
    }
  })
}
print(table(find_min_nonzero_indices(loglikehood_total1)))