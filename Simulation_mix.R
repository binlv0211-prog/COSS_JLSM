library(pgdraw)
library(LaplacesDemon)
library(MASS)
library(pROC)
library(foreach)
library(doParallel)
source("C:/Users/33701/Desktop/code/network_only.R")
source("C:/Users/33701/Desktop/code/network_mix.R")
source("C:/Users/33701/Desktop/code/Y_only_mix.R")
source("C:/Users/33701/Desktop/code/data_process.R")


nrun = 15000
burn = 10000
thin = 5
cn = c(50,100,150,300)
cq = c(5,10,15,30)
k = 3
alpha_H = 8
a_sig = 1
b_sig = 1
a_theta = 2
b_theta = 2
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
b1_l = -3
b1_u = -1.5
b2_l = -1.5
b2_u = 0
b3_l = 0
b3_u = 1.5
b4_l = 1.5
b4_u = 3
pll = 50
replation = 100

cl = makeCluster(pll)      
registerDoParallel(cl)
res = foreach(iter = 1:replation, .verbose = TRUE, .packages = c("MASS","pgdraw","LaplacesDemon","truncnorm"), .combine = list, .multicombine = TRUE) %dopar% {
  out_net = list()
  out_net_Y = list()
  out_Y = list()
  data_ture = list()  
  for (i in 1:length(cn)) {
    n = cn[i]
    q = cq[i]
    
    delta_n = 1 + 0.001 * n 
    alpha_ture = get_alpha(n,alpha_l,alpha_u)
    Z_ture = get_Z(n,k)
    A = get_A(alpha_ture,Z_ture)
    gamma1_ture = get_gamma(q)
    B1_ture = get_B(q,k,B_l,B_u,T)
    Y1 = get_Y(gamma1_ture,Z_ture,B1_ture,continous = T)
    gamma2_ture = get_gamma(q)
    B2_ture = get_B(q,k,B_l,B_u,T)
    Y2 = get_Y(gamma2_ture,Z_ture,B2_ture,continous = F)
    b1<-runif(q,-3,-1.5)
    b2<-runif(q,-1.5,0)
    b3<-runif(q,0,1.5)
    b4<-runif(q,1.5,3)
    B3_ture = get_B(q,k,B_l,B_u,T)
    Y3 = get_item_Y(Z_ture,B3_ture,b1,b2,b3,b4)   
    data_ture[[i]] = list("alpha" = alpha_ture,"Z" = Z_ture,"A" = A,"gamma1" = gamma1_ture,"B1" = B1_ture,"Y1" = Y1,
                          "gamma2" = gamma2_ture,"B2" = B2_ture,"Y2" = Y2,
                          "b1" = b1,"b2" = b2,"b3" = b3,"b4" = b4,"B3" = B3_ture,"Y3" = Y3)
    net_out = network_only_1123(A, 1314, nrun,burn,thin, delta_n, alpha_H,a_theta, b_theta, 
                                theta_inf, start_adapt, Hmax, alpha0,alpha1)
    net_Y_out = network_mix(A,Y1,Y2,Y3,nrun,burn,thin,delta_n,alpha_H,a_theta,b_theta,a_sig,b_sig,
                            theta_inf,Hmax,start_adapt,alpha0,alpha1)    
    Y_out = Y_mix(Y1,Y2,Y3,nrun,burn,thin,delta_n,alpha_H,a_sig,b_sig,a_theta_B,b_theta_B,
                  theta_inf = 0.05,start_adapt,Hmax,alpha0,alpha1)
    
    out_net[[i]] = net_out
    out_net_Y[[i]] = net_Y_out
    out_Y[[i]] = Y_out
  }
  out = list("net" = out_net,"net_Y" = out_net_Y,"Y" = out_Y)
  re = list("out" = out,"ture" = data_ture)
  return(re)  
}
stopCluster(cl)