library(pgdraw)
library(LaplacesDemon)
library(MASS)
library(pROC)
library(foreach)
library(doParallel)
source("C:/Users/33701/Desktop/code/network_only.R")
source("C:/Users/33701/Desktop/code/network_binary.R")
source("C:/Users/33701/Desktop/code/network_normal.R")
source("C:/Users/33701/Desktop/code/Y_only_binary.R")
source("C:/Users/33701/Desktop/code/Y_only_normal.R")
source("C:/Users/33701/Desktop/code/network_normal_kl.R")
source("C:/Users/33701/Desktop/code/network_binary_kl.R")
source("C:/Users/33701/Desktop/code/data_process.R")

nrun = 15000
burn = 10000
thin = 5
cn = c(50,100,150,300)
cq = c(10,20,30,60)
k = 3
alpha_H_c = c(2.5,5,10,5,5,5)
a_sig = 1
b_sig = 1
a_theta_c = c(3,3,3,1.5,3,1.5)
b_theta_c = c(3,3,3,3,1.5,1.5)
a_theta_B = 2
b_theta_B = 2
theta_inf = 0.10
start_adapt = 500
Hmax = 7
alpha0 = 1
alpha1 = 5*10^(-4)
alpha_l = -0.5
alpha_u = 0.5#保证网络的密度大于0.5
B_l = 0.25
B_u = 1.25
pll = 50
replation = 100
ni = 1
n = cn[ni]
q = cq[ni]

library(foreach)
library(doParallel)
cl = makeCluster(pll)      
registerDoParallel(cl)
res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon"), .combine = list, .multicombine = TRUE) %dopar% {
  out_net_Y = list()
  data_ture = list()
  for (i in 1:length(alpha_H_c)) {
    alpha_H = alpha_H_c[i]
    a_theta = a_theta_c[i]
    b_theta = b_theta_c[i]
    delta_n = 1 + 0.001 * n 
    alpha_ture = get_alpha(n,alpha_l,alpha_u)
    Z_ture = get_Z(n,k)
    A = get_A(alpha_ture,Z_ture)
    gamma_ture = get_gamma(q)
    B_ture = get_B(q,k,B_l,B_u,T)
    Y = get_Y(gamma_ture,Z_ture,B_ture,continous = T)
    data_ture[[i]] = list("alpha" = alpha_ture,"Z" = Z_ture,"gamma" = gamma_ture,"B" = B_ture,"A" = A,"Y" = Y)
    net_Y_out = network_nomal(A, Y,nrun,burn,thin, delta_n, alpha_H,a_sig,b_sig
                              ,a_theta,b_theta, #a_theta_B,b_theta_B,
                              theta_inf,start_adapt, Hmax, alpha0,alpha1)
    out_net_Y[[i]] = net_Y_out
  }
  out = list("net_Y" = out_net_Y)
  re = list("out" = out,"ture" = data_ture)
  return(re)
}
stopCluster(cl)


cl = makeCluster(pll)      
registerDoParallel(cl)
res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon"), .combine = list, .multicombine = TRUE) %dopar% {
  out_net_Y = list()
  data_ture = list()
  for (i in 1:length(alpha_H_c)) {
    alpha_H = alpha_H_c[i]
    a_theta = a_theta_c[i]
    b_theta = b_theta_c[i]
    delta_n = 1 + 0.001 * n 
    alpha_ture = get_alpha(n,alpha_l,alpha_u)
    Z_ture = get_Z(n,k)
    A = get_A(alpha_ture,Z_ture)
    gamma_ture = get_gamma(q)
    B_ture = get_B(q,k,B_l,B_u,T)
    Y = get_Y(gamma_ture,Z_ture,B_ture,continous = F)
    data_ture[[i]] = list("alpha" = alpha_ture,"Z" = Z_ture,"gamma" = gamma_ture,"B" = B_ture,"A" = A,"Y" = Y)
    net_Y_out = network_briny_1130(A,Y,nrun,burn,thin,delta_n,alpha_H,a_theta,b_theta,
                                   theta_inf,start_adapt,Hmax,alpha0,alpha1)
    out_net_Y[[i]] = net_Y_out
  }
  out = list("net_Y" = out_net_Y)
  re = list("out" = out,"ture" = data_ture)
  return(re)
}
stopCluster(cl)