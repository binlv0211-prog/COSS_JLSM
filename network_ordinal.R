library(pgdraw)
library(MASS)
library(truncnorm)
library(LaplacesDemon)
joint_item_Y = function(A,Y,nrun,burn,thin,delta_n,alpha_H,a_theta,b_theta,
                        theta_inf,Hmax,start_adapt,alpha0,alpha1){
  n = dim(Y)[1]
  q = dim(Y)[2]
  u<-runif(nrun)
  #初始化
  H = Hmax + 1
  Hstar = Hmax
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
  theta_inv = rep(1,H)
  theta_inv_B = rep(1,H)
  w = rep(1,H)
  zta = rep(1,H)
  B = matrix(rnorm(H * q),nrow = H, ncol = q)
  for (i in 2:H) {
    for (j in 1:(i-1)) {
      B[i,j] = 0
    }
  }
  Phi = matrix(1,n,q)#文章中的Z
  W = matrix(1,n,q)
  b1<-rtruncnorm(q, a=-3, b=-1.5, mean = 0, sd = 1)
  b2<-rtruncnorm(q, a=-1.5, b=0, mean = 0, sd = 1)
  b3<-rtruncnorm(q, a=0, b=1.5, mean = 0, sd = 1)
  b4<-rtruncnorm(q, a=1.5, b=3, mean = 0, sd = 1)
  #存数据
  N_sample = ceiling((nrun - burn)/thin)
  H_hat = rep(NA,nrun)
  Z_hat = list()
  B_hat = list()
  alpha_hat = matrix(0,N_sample,n)
  b1_hat = matrix(0,N_sample,q)
  b2_hat = matrix(0,N_sample,q)
  b3_hat = matrix(0,N_sample,q)
  b4_hat = matrix(0,N_sample,q)
  m = 1
  for (run in 1:nrun) {
    #sample Phi
    Theta = Z %*% B
    x0 = matrix(-10,n,q)
    b1_temp = outer(rep(1,n),b1)
    b2_temp = outer(rep(1,n),b2)
    b3_temp = outer(rep(1,n),b3)
    b4_temp = outer(rep(1,n),b4)
    x1 = matrix(10,n,q)
    phi1<-matrix(rtruncnorm(1, a=x0, b=b1_temp, mean = Theta, sd = sqrt(1/W)),n,q)#1
    phi1 = ifelse(Y==1,phi1,0)
    phi2<-matrix(rtruncnorm(1, a=b1_temp, b=b2_temp, mean = Theta, sd = sqrt(1/W)),n,q)#2
    phi2 = ifelse(Y==2,phi2,0)
    phi3<-matrix(rtruncnorm(1, a=b2_temp, b=b3_temp, mean = Theta, sd = sqrt(1/W)),n,q)#3
    phi3 = ifelse(Y==3,phi3,0)
    phi4<-matrix(rtruncnorm(1, a=b3_temp, b=b4_temp, mean = Theta, sd = sqrt(1/W)),n,q)#4
    phi4 = ifelse(Y==4,phi4,0)
    phi5<-matrix(rtruncnorm(1, a=b4_temp, b=x1, mean = Theta, sd = sqrt(1/W)),n,q)#5
    phi5 = ifelse(Y==5,phi5,0)
    Phi = phi1 + phi2 + phi3 + phi4 + phi5
    #sample W
    Theta_w = Z %*% B - Phi
    W = matrix(pgdraw(2,Theta_w),n,q) 
    #sample B
    for (j in 1:q) {
      W_j = W[,j]
      phi_j = Phi[,j]
      if(j<H){
        sigma_Bj = chol2inv(chol(diag(theta_inv_B[1:j],nrow = j)+ t(Z[,1:j]) %*% diag(W_j,nrow = n) %*% Z[,1:j]))
        u_Bj = sigma_Bj %*% t(Z[,1:j]) %*%diag(W_j,nrow = n)%*%phi_j
        B[1:j,j] = mvrnorm(1, u_Bj, sigma_Bj)
      }
      else{
        sigma_Bj = chol2inv(chol(t(Z)%*%diag(W_j,nrow = n)%*%Z + diag(theta_inv_B,nrow = H)))
        u_Bj = sigma_Bj%*%t(Z)%*%diag(W_j,nrow = n)%*%phi_j
        B[,j] = mvrnorm(1,u_Bj,sigma_Bj)
      }
    }
    # for (j in 1:q) {
    #   W_j = W[,j]
    #   phi_j = Phi[,j]
    #   sigma_Bj = chol2inv(chol(t(Z)%*%diag(W_j,nrow = n)%*%Z + diag(H)))
    #   u_Bj = sigma_Bj%*%t(Z)%*%diag(W_j,nrow = n)%*%phi_j
    #   B[,j] = mvrnorm(1,u_Bj,sigma_Bj)
    # }
    #sample D
    theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
    D_A = matrix(pgdraw(1, theta_A), nrow = n,ncol = n)
    D_A = (D_A + t(D_A)) / 2
    #update alpha
    Z_temp = Z %*% t(Z)
    for(i in 1:n){
      sigma_alphai = 1 / (sum(D_A[i,]) - D_A[i,i] + 1/100)#此处默认alpha先验的方差为1，后续可能会改动
      u_temp = A[i,] - 0.5 - D_A[i,] * (alpha + Z_temp[i,])
      u_alphai = sigma_alphai * (sum(u_temp) - u_temp[i])
      alpha[i] = rnorm(1,u_alphai,sqrt(sigma_alphai))
    }
    for (i in 1:n) {
      W_i = W[i,]
      phi_i = Phi[i,]
      Z_i = Z[-i,]
      D_Ai = (D_A[i,])[-i]
      alp_cons = alpha[i] + alpha[-i]
      kappa_Ai = (A[i,])[-i] - 0.5
      sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i + B%*%diag(W_i,nrow = q)%*%t(B)))
      u_Zi = sigma_Zi%*%((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) + (B%*%diag(W_i,nrow = q)%*%phi_i))
      Z[i,] = mvrnorm(1,u_Zi,sigma_Zi)
    }
    # sample zta
    lhd_spike<-rep(0,H)
    lhd_slab<-rep(0,H)
    for(h in 1:H){
      lhd_spike[h] = exp(sum(log(dnorm(Z[,h], mean = 0, sd = theta_inf^(1/2), log = FALSE))))
      #lhd_spike[h] = min(G, lhd_spike[h])
      lhd_slab[h] = dmvt(x = Z[,h], mu=rep(0,n), S=(b_theta/a_theta)*diag(n), df=2*a_theta)
      prob_h = w*c(rep(lhd_spike[h],h),rep(lhd_slab[h],H - h))
      # print(lhd_spike[h])
      # print(lhd_slab[h])
      # print("----")
      if (sum(prob_h) == 0){
        prob_h = c(rep(0,H-1),1)
      }
      else{
        prob_h = prob_h/sum(prob_h)
      }
      zta[h] = c(1:H)%*%rmultinom(n=1, size=1, prob=prob_h)
    }
    #sample v and update w
    v = rep(NA,H)
    for (h in 1:(H - 1)){
      if (h == 1){
        v[h] = rbeta(1, shape1 = Hmax**(delta_n) + sum(zta == h), shape2 = 1 + sum(zta > h))
      }else{
        v[h] = rbeta(1, shape1 = alpha_H + sum(zta == h), shape2 = 1 + sum(zta > h))
      }
    }
    v[H] = 1
    w[1] = v[1]
    for (h in 2:H){
      w[h] = v[h]*prod(1-v[1:(h-1)])  
    }
    # 6) sample theta^{-1}
    for (h in 1:H){
      if (zta[h] <= h){
        theta_inv[h] = theta_inf^(-1)
      }
      else{
        theta_inv[h] = rgamma(n=1,shape = a_theta + 0.5 * n,rate=b_theta + 0.5 * t(Z[,h]) %*% Z[,h])
      }
    }
    #sample b1234
    pb1<-ifelse(Y==1,Phi,NA)
    pb2<-ifelse(Y==2,Phi,NA)
    pb3<-ifelse(Y==3,Phi,NA)
    pb4<-ifelse(Y==4,Phi,NA)
    pb5<-ifelse(Y==5,Phi,NA)
    pb11<-apply(pb1,2,max,na.rm=TRUE)
    pb20<-apply(pb2,2,min,na.rm=TRUE)
    pb21<-apply(pb2,2,max,na.rm=TRUE)
    pb30<-apply(pb3,2,min,na.rm=TRUE)
    pb31<-apply(pb3,2,max,na.rm=TRUE)
    pb40<-apply(pb4,2,min,na.rm=TRUE)
    pb41<-apply(pb4,2,max,na.rm=TRUE)
    pb50<-apply(pb5,2,min,na.rm=TRUE)
    pb11<-pmax(pb11,-3)
    pb20<-pmin(pb20,-1.5)
    pb21<-pmax(pb21,-1.5)
    pb30<-pmin(pb30,0)
    pb31<-pmax(pb31,0)
    pb40<-pmin(pb40,1.5)
    pb41<-pmax(pb41,1.5)
    pb50<-pmin(pb50,3)
    Q<-matrix(c(pb11,pb20),ncol=2)
    pb11 = apply(Q, 1, min)
    pb20 = apply(Q, 1, max)
    Q<-matrix(c(pb21,pb30),ncol=2)
    pb21 = apply(Q, 1, min)
    pb30 = apply(Q, 1, max)
    Q<-matrix(c(pb31,pb40),ncol=2)
    pb31 = apply(Q, 1, min)
    pb40 = apply(Q, 1, max)
    Q<-matrix(c(pb41,pb50),ncol=2)
    pb41 = apply(Q, 1, min)
    pb50 = apply(Q, 1, max)
    b10 = runif(q,pb11,pb20)
    b20 = runif(q,pb21,pb30)
    b30 = runif(q,pb31,pb40)
    b40 = runif(q,pb41,pb50)
    b1 = rtruncnorm(1, a=pb11, b=pb20, mean = b10, sd = 10^5)
    b2 = rtruncnorm(1, a=pb21, b=pb30, mean = b20, sd = 10^5)
    b3 = rtruncnorm(1, a=pb31, b=pb40, mean = b30, sd = 10^5)
    b4 = rtruncnorm(1, a=pb41, b=pb50, mean = b40, sd = 10^5)
    #update H[t]
    active = which(zta > c(1:H))
    Hstar = length(active)
    if (run >= start_adapt & u[run] <= exp(alpha0 + alpha1 * run)){
      if (Hstar < H - 1){
        # set truncation to Hstar[t] and subset all variables, keeping only active columns
        H = Hstar + 1
        theta_inv = c(theta_inv[active],theta_inf^(-1))
        w = c(w[active],1-sum(w[active]))
        Z = cbind(Z[,active],rnorm(n,mean=0,sd=sqrt(theta_inf)))
        B = rbind(B[active,],c(rep(0,H-1),rnorm(q-H+1)))
        #theta_inv_B = c(theta_inv_B[active],rgamma(1,a_theta_B,b_theta_B))
        #B = rbind(B[active,],rnorm(q,0,sqrt(theta_inv_B[H])))
        zta = c(zta[active],H-1)
      } else if (H < Hmax) {
        # increase truncation by 1 and extend all variables, sampling from the prior/model
        H = H + 1
        v[H - 1] = rbeta(1,shape1=alpha_H,shape2=1)
        v = c(v,1)
        w = rep(NA,H)
        w[1] = v[1]
        for (h in 2:H){
          w[h] = v[h]*prod(1-v[1:(h-1)])
        }
        theta_inv = c(theta_inv,theta_inf^(-1))
        Z = cbind(Z,rnorm(n,mean=0,sd=sqrt(theta_inf)))
        B = rbind(B,c(rep(0,H-1),rnorm(q-H+1)))
        #theta_inv_B = c(theta_inv_B,rgamma(1,a_theta_B,b_theta_B))
        #B = rbind(B,rnorm(q,0,sqrt(theta_inv_B[H])))
        zta = c(zta,H-1)
      }
    }
    H_hat[run] = Hstar
    if((run > burn) &((run-burn) %% thin == 0)){
      b1_hat[m,] = b1
      b2_hat[m,] = b2
      b3_hat[m,] = b3
      b4_hat[m,] = b4
      alpha_hat[m,] = alpha
      if(Hstar>0){
        B_hat[[m]] = B[1:Hstar,, drop=FALSE]
        Z_hat[[m]] = Z[,1:Hstar, drop=FALSE]
      }else{
        B_hat[[m]] = B
        Z_hat[[m]] = Z
      }
      m = m + 1
    }    
  }
  output = list("H" = H_hat,"Z" = Z_hat,"B" = B_hat,"alpha" = alpha_hat,"b1" = b1_hat,
                "b2" = b2_hat,"b3" = b3_hat,"b4" = b4_hat)
  return(output)    
}
# alpha_H = 5
# a_sig = 1
# b_sig = 1
# a_theta = 2
# b_theta = 2
# a_theta_B = 2
# b_theta_B = 2
# theta_inf = 0.05
# start_adapt = 100
# Hmax = 10
# alpha0 = 1
# alpha1 = 5*10^(-4)
# alpha_l = -0.5
# alpha_u = 0.5#保证网络的密度大于0.5
# B_l = 1
# B_u = 2
# library(truncnorm)
# source("/home/binlv/data_process.R")
# n = 400
# H = 6
# q = 20
# alpha = get_alpha(n,alpha_l,alpha_u)
# Z = matrix(rnorm(n * H),n,H)
# A = get_A(alpha,Z)
# B = get_B(q,H,1,2,T)
# #B = matrix(rnorm(q * H),H,q)
# b1<-rtruncnorm(q, a=-3, b=-1.5, mean = 0, sd = 1)
# b2<-rtruncnorm(q, a=-1.5, b=0, mean = 0, sd = 1)
# b3<-rtruncnorm(q, a=0, b=1.5, mean = 0, sd = 1)
# b4<-rtruncnorm(q, a=1.5, b=3, mean = 0, sd = 1)
# b1_temp = outer(rep(1,n),b1)
# b2_temp = outer(rep(1,n),b2)
# b3_temp = outer(rep(1,n),b3)
# b4_temp = outer(rep(1,n),b4)  
# Psi = array(0,c(n,q,6))
# P = array(0,c(n,q,5))
# Psi[,,1] = 1
# Psi[,,2] = exp(Z%*%B - b1_temp)/(1 + exp(Z%*%B - b1_temp))
# Psi[,,3] = exp(Z%*%B - b2_temp)/(1 + exp(Z%*%B - b2_temp))
# Psi[,,4] = exp(Z%*%B - b3_temp)/(1 + exp(Z%*%B - b3_temp))
# Psi[,,5] = exp(Z%*%B - b4_temp)/(1 + exp(Z%*%B - b4_temp))
# P[,,1] = Psi[,,1] - Psi[,,2]
# P[,,2] = Psi[,,2] - Psi[,,3]
# P[,,3] = Psi[,,3] - Psi[,,4]
# P[,,4] = Psi[,,4] - Psi[,,5]
# P[,,5] = Psi[,,5] - Psi[,,6]
# Y = matrix(0,n,q)
# for (i in 1:n) {
#   for (j in 1:q) {
#     Y[i,j] = sample(1:5,1,replace = T,prob = P[i,j,])
#   }
# }
# out = joint_item_Y(A,Y,15000,10000,5,alpha_H,a_theta,b_theta,
#                    theta_inf,Hmax,start_adapt,alpha0,alpha1)