library(pgdraw)
library(LaplacesDemon)
library(MASS)
library(foreach)
library(doParallel)
source("functions/real_data3_functions.R")
# 为每个数据文件创建单独的环境
friendship_env <- new.env()
demographic_env <- new.env()
substances_env <- new.env()
various_env <- new.env()
lifestyle_env <- new.env()
geographic_env <- new.env()
selections_env <- new.env()
# 加载到各自的环境
load("datas/Glasgow-friendship.RData", envir = friendship_env)
load("datas/Glasgow-demographic.RData", envir = demographic_env)
load("datas/Glasgow-substances.RData", envir = substances_env)
load("datas/Glasgow-various.RData", envir = various_env)
load("datas/Glasgow-lifestyle.RData", envir = lifestyle_env)
load("datas/Glasgow-geographic.RData", envir = geographic_env)
load("datas/Glasgow-selections.RData", envir = selections_env)
selections = selections_env$selection129
n = sum(selections)

F1 = friendship_env$friendship.1[selections,selections]
A = matrix(0,n,n)
for (i in 1:n) {
  for (j in 1:n) {
    if(F1[i,j] + F1[j,i] > 0){A[i,j] = 1}
  }
}
Y1 = matrix(substances_env$alcohol[selections,1],n,1) # k=5
q1 = dim(Y1)[2]
MissY1 = matrix(0,n,q1)
MissY1[is.na(Y1)] = 1
Y1[is.na(Y1)] = 1
Y2 = matrix(substances_env$tobacco[selections,1],n,1)
Y3 = cbind(substances_env$cannabis[selections,1],lifestyle_env$leisure1[selections,]) # k=4
q3 = dim(Y3)[2]
MissY3 = matrix(0,n,q3)
MissY3[is.na(Y3)] = 1
Y3[is.na(Y3)] = 1
Y4 = lifestyle_env$music1[selections,]


pll = 50
replation = 100
cl = makeCluster(pll)      
registerDoParallel(cl)
res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon","truncnorm"), .combine = list, .multicombine = TRUE) %dopar% {
  net_Y_out = network_mix_real_data(A,Y1,Y2,Y3,Y4,MissY1,MissY3,nrun = 15000,burn = 10000,
                                    thin = 5, delta_n = 0.001 * n, alpha_H = 3,a_theta = 1,b_theta = 1
                                    ,a_sig = 3,b_sig = 3, theta_inf = 0.1,Hmax = 10,start_adapt = 500
                                    ,alpha0 = -1, alpha1 = -5e-4)
  return(net_Y_out)
}
stopCluster(cl)


cl = makeCluster(pll)      
registerDoParallel(cl)
res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon","truncnorm"), .combine = list, .multicombine = TRUE) %dopar% {
  res_criterion = list()
  res_criterion[[1]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,1,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[2]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,2,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[3]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,3,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[4]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,4,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[5]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,5,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[6]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,6,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[7]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,7,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[8]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,8,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[9]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,9,nrun = 15000,burn = 10000,thin = 5)
  res_criterion[[10]] = network_mix_real_data_kf(A,Y1,Y2,Y3,Y4,MissY1,MissY3,10,nrun = 15000,burn = 10000,thin = 5)
  return(res_criterion)
}
stopCluster(cl)

K = 5
cl = makeCluster(pll)      
registerDoParallel(cl)
res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon","truncnorm"), .combine = list, .multicombine = TRUE) %dopar% {
  res_criterion = list()
  loglikehood = matrix(0,10,K)
  folds <- create_folds(n, k = K)  
  for(i in 1:K){
    train_ids <- folds[[i]]$train
    test_ids <- folds[[i]]$test
    n_train = length(train_ids)
    n_test = length(test_ids)
    A_train = A[train_ids,train_ids,drop=FALSE]
    Y1_train = Y1[train_ids,,drop=FALSE]
    MissY1_train = MissY1[train_ids,,drop=FALSE]
    Y2_train = Y2[train_ids,,drop=FALSE]
    Y3_train = Y3[train_ids,,drop=FALSE]
    MissY3_train = MissY3[train_ids,,drop=FALSE]
    Y4_train = Y4[train_ids,,drop=FALSE]
    A_test = A[test_ids,test_ids,drop=FALSE]
    Y1_test = Y1[test_ids,,drop=FALSE]
    MissY1_test = MissY1[test_ids,,drop=FALSE]
    Y2_test = Y2[test_ids,,drop=FALSE]
    Y3_test = Y3[test_ids,,drop=FALSE]
    MissY3_test = MissY3[test_ids,,drop=FALSE]
    Y4_test = Y4[test_ids,,drop=FALSE]
    
    for (h in 1:10) {
      kf_out = network_mix_real_data_kf(A_train,Y1_train,Y2_train,Y3_train,Y4_train,MissY1_train,MissY3_train,h,nrun = 15000,burn = 10000,thin = 5)
      b11 = apply(kf_out$b11, 2, mean)
      b12 = apply(kf_out$b12, 2, mean)
      b13 = apply(kf_out$b13, 2, mean)
      b14 = apply(kf_out$b14, 2, mean)
      b21 = apply(kf_out$b21, 2, mean)
      b22 = apply(kf_out$b22, 2, mean)
      b31 = apply(kf_out$b31, 2, mean)
      b32 = apply(kf_out$b32, 2, mean)
      b33 = apply(kf_out$b33, 2, mean)
      gamma4 = apply(kf_out$gamma4, 2, mean)
      B1 = apply(kf_out$B1,c(2,3),mean)
      B2 = apply(kf_out$B2,c(2,3),mean)
      B3 = apply(kf_out$B3,c(2,3),mean)
      B4 = apply(kf_out$B4,c(2,3),mean) 
      
      kf_outz = network_mix_real_data_kf_getZ(A_test,Y1_test,Y2_test,Y3_test,Y4_test,MissY1_test,MissY3_test,B1,B2,B3,B4,b11,b12,b13,b14,b21,b22,b31,b32,b33,gamma4,h,nrun = 15000,burn = 10000,thin = 5)
      alpha = apply(kf_outz$alpha,2,mean)
      Z = apply(kf_outz$Z,c(2,3),mean)      
      loglikehood[h,i] = get_loglikehood_data3(A_test,Y1_test,Y2_test,Y3_test,Y4_test,MissY1_test,MissY3_test,b11,b12,b13,b14,b21,
                                               b22,b31,b32,b33,gamma4,B1,B2,B3,B4,alpha,Z)
    }
  }
  return(loglikehood)
}
stopCluster(cl)