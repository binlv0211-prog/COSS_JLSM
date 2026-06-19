library(pgdraw)
library(LaplacesDemon)
library(MASS)
library(pROC)
library(foreach)
library(doParallel)
source("functions/network_only.R")
source("functions/network_binary.R")
source("functions/network_normal.R")
source("functions/Y_only_binary.R")
source("functions/Y_only_normal.R")
source("functions/network_normal_kl.R")
source("functions/network_binary_kl.R")
source("functions/data_process.R")

nrun = 15000
burn = 10000
thin = 5
n = 100
q = 20
k = 3
a_sig = 1
b_sig = 1
a_theta = 3
b_theta = 3
a_theta_B = 3
b_theta_B = 3
theta_inf = 0.10
start_adapt = 500
Hmax = 7
alpha_H = 8
alpha0 = 1
alpha1 = 5*10^(-4)
alphalc = seq(from = -3, to = -0.375, length.out = 8)
alphauc = seq(from = -1, to = -0.125, length.out = 8) # 0.12 - 0.42

B_l = 0.25
B_u = 1.25
pll = 1
replation = 100


res_normal = list()
for (i in 1:length(alphauc)){

  alpha_l = alphalc[i]
  alpha_u = alphauc[i]#   
  
  cl = makeCluster(pll)      
  registerDoParallel(cl)
  
  res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon"), .combine = list, .multicombine = TRUE) %dopar% {
    alpha_l = alphalc[i]
    alpha_u = alphauc[i]#
    out = list()
    delta_n = 1 + 0.001 * n 
    alpha_ture = get_alpha(n,alpha_l,alpha_u)
    Z_ture = get_Z(n,k)
    A = get_A(alpha_ture,Z_ture)
    gamma_ture = get_gamma(q)
    B_ture = get_B(q,k,B_l,B_u,T)  
    Y = get_Y(gamma_ture,Z_ture,B_ture,continous = T)
    data_ture = list("alpha" = alpha_ture,"Z" = Z_ture,"gamma" = gamma_ture,"B" = B_ture,"A" = A,"Y" = Y)
    net_out = network_only_1123(A, 1314, nrun,burn,thin, delta_n, alpha_H,a_theta, b_theta, 
                                theta_inf, start_adapt, Hmax, alpha0,alpha1)
    net_Y_out = network_nomal(A, Y,nrun,burn,thin, delta_n,alpha_H,a_sig,b_sig
                              ,a_theta,b_theta, #a_theta_B,b_theta_B,
                              theta_inf,start_adapt, Hmax,alpha0,alpha1)
    out = list("net" = net_out,"net_Y" = net_Y_out)
    re = list("true" = B_ture, "out" = out)
    return(re)
  }
  stopCluster(cl)
  res_normal[[i]] = res
}

res_binary = list()
for (i in 1:length(alphauc)){
  
  alpha_l = alphalc[i]
  alpha_u = alphauc[i]#   

  cl = makeCluster(pll)      
  registerDoParallel(cl)
  
  res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon"), .combine = list, .multicombine = TRUE) %dopar% {
    out = list()
    delta_n = 1 + 0.001 * n 
    alpha_ture = get_alpha(n,alpha_l,alpha_u)
    Z_ture = get_Z(n,k)
    A = get_A(alpha_ture,Z_ture)
    gamma_ture = get_gamma(q)
    B_ture = get_B(q,k,B_l,B_u,T)  
    Y = get_Y(gamma_ture,Z_ture,B_ture,continous = F)
    data_ture = list("alpha" = alpha_ture,"Z" = Z_ture,"gamma" = gamma_ture,"B" = B_ture,"A" = A,"Y" = Y)
    net_out = network_only_1123(A, 1313, nrun,burn,thin, delta_n, alpha_H,a_theta, b_theta, 
                                theta_inf, start_adapt, Hmax,  alpha0,alpha1)
    net_Y_out = network_briny_1130(A,Y,nrun,burn,thin,delta_n,alpha_H,a_theta,b_theta,
                                   theta_inf,start_adapt,Hmax,alpha0,alpha1)
    out = list("net" = net_out,"net_Y" = net_Y_out)
    re = list("true" = B_ture, "out" = out)
    return(re)
  }
  stopCluster(cl)
  res_binary[[i]] = res
}