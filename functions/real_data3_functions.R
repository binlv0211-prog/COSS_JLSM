library(pgdraw)
library(MASS)
library(truncnorm)
library(LaplacesDemon)

safe_max <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(-10) #
  return(max(x_clean))
}

safe_min <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(10) # 
  return(min(x_clean))
}


get_item_Y1 = function(Z,B,b1,b2,b3,b4){
  n = dim(Z)[1]
  q = dim(B)[2]
  b1_temp = outer(rep(1,n),b1)
  b2_temp = outer(rep(1,n),b2)
  b3_temp = outer(rep(1,n),b3)
  b4_temp = outer(rep(1,n),b4)  
  Psi = array(0,c(n,q,6))
  P = array(0,c(n,q,5))
  Psi[,,1] = 1
  Psi[,,2] = exp(Z%*%B - b1_temp)/(1 + exp(Z%*%B - b1_temp))
  Psi[,,3] = exp(Z%*%B - b2_temp)/(1 + exp(Z%*%B - b2_temp))
  Psi[,,4] = exp(Z%*%B - b3_temp)/(1 + exp(Z%*%B - b3_temp))
  Psi[,,5] = exp(Z%*%B - b4_temp)/(1 + exp(Z%*%B - b4_temp))
  P[,,1] = Psi[,,1] - Psi[,,2]
  P[,,2] = Psi[,,2] - Psi[,,3]
  P[,,3] = Psi[,,3] - Psi[,,4]
  P[,,4] = Psi[,,4] - Psi[,,5]
  P[,,5] = Psi[,,5] - Psi[,,6]
  Y = matrix(0,n,q)
  for (i in 1:n) {
    for (j in 1:q) {
      Y[i,j] = sample(1:5,1,replace = T,prob = P[i,j,])
    }
  }
  return(Y)
}

get_item_Y3 = function(Z,B,b1,b2,b3){
  n = dim(Z)[1]
  q = dim(B)[2]
  b1_temp = outer(rep(1,n),b1)
  b2_temp = outer(rep(1,n),b2)
  b3_temp = outer(rep(1,n),b3)
  Psi = array(0,c(n,q,5))
  P = array(0,c(n,q,4))
  Psi[,,1] = 1
  Psi[,,2] = exp(Z%*%B - b1_temp)/(1 + exp(Z%*%B - b1_temp))
  Psi[,,3] = exp(Z%*%B - b2_temp)/(1 + exp(Z%*%B - b2_temp))
  Psi[,,4] = exp(Z%*%B - b3_temp)/(1 + exp(Z%*%B - b3_temp))
  P[,,1] = Psi[,,1] - Psi[,,2]
  P[,,2] = Psi[,,2] - Psi[,,3]
  P[,,3] = Psi[,,3] - Psi[,,4]
  P[,,4] = Psi[,,4] - Psi[,,5]
  Y = matrix(0,n,q)
  for (i in 1:n) {
    for (j in 1:q) {
      Y[i,j] = sample(1:4,1,replace = T,prob = P[i,j,])
    }
  }
  return(Y)
}
# Y1 is five item
# Y2 is  three item with missing
# Y3 is four item with missing
# Y4 is binary

network_mix_real_data = function(A,Y1,Y2,Y3,Y4,MissY1,MissY3,nrun,burn,thin,delta_n,alpha_H,a_theta,b_theta,a_sig,b_sig,
                                 theta_inf,Hmax,start_adapt,alpha0,alpha1){
  n = dim(Y1)[1]
  q1 = dim(Y1)[2]
  q2 = dim(Y2)[2]
  q3 = dim(Y3)[2]
  q4 = dim(Y4)[2]
  u<-runif(nrun)
  H = Hmax + 1
  Hstar = Hmax
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)  
  theta_inv = rep(1,H)  
  
  theta_inv_B1 = rep(1,H)
  theta_inv_B2 = rep(1,H)
  theta_inv_B3 = rep(1,H)
  theta_inv_B4 = rep(1,H)
  w = rep(1,H)
  zta = rep(1,H)    
  B1 = matrix(rnorm(H * q1),nrow = H, ncol = q1)
  B2 = matrix(rnorm(H * q2),nrow = H, ncol = q2)
  B3 = matrix(rnorm(H * q3),nrow = H, ncol = q3)
  B4 = matrix(rnorm(H * q4),nrow = H, ncol = q4)
  
  Phi1 = matrix(1,n,q1)#
  W1 = matrix(1,n,q1)
  b11<-rtruncnorm(q1, a=-3, b=-1.5, mean = 0, sd = 1)
  b12<-rtruncnorm(q1, a=-1.5, b=0, mean = 0, sd = 1)
  b13<-rtruncnorm(q1, a=0, b=1.5, mean = 0, sd = 1)
  b14<-rtruncnorm(q1, a=1.5, b=3, mean = 0, sd = 1)  
  
  Phi2 = matrix(1,n,q2)#
  W2 = matrix(1,n,q2)
  b21<-rtruncnorm(q2, a=-3, b=0, mean = 0, sd = 1)
  b22<-rtruncnorm(q2, a=0, b=3, mean = 0, sd = 1)
  
  Phi3 = matrix(1,n,q3)#
  W3 = matrix(1,n,q3)
  b31<-rtruncnorm(q3, a=-3, b=-1, mean = 0, sd = 1)
  b32<-rtruncnorm(q3, a=-1, b=1, mean = 0, sd = 1)
  b33<-rtruncnorm(q3, a=1, b=3, mean = 0, sd = 1)
  
  gamma4 = rep(1,q4)
  
  N_sample = ceiling((nrun - burn)/thin)
  H_hat = rep(NA,nrun)  
  
  Z_hat = list()
  B1_hat = list()
  B2_hat = list()
  B3_hat = list()
  B4_hat = list()
  
  alpha_hat = matrix(0,N_sample,n)  
  
  b11_hat = matrix(0,N_sample,q1)
  b12_hat = matrix(0,N_sample,q1)
  b13_hat = matrix(0,N_sample,q1)
  b14_hat = matrix(0,N_sample,q1)
  
  b21_hat = matrix(0,N_sample,q2)
  b22_hat = matrix(0,N_sample,q2)
  
  b31_hat = matrix(0,N_sample,q3)
  b32_hat = matrix(0,N_sample,q3)
  b33_hat = matrix(0,N_sample,q3)
  
  gamma4_hat = matrix(0,N_sample,q4)
  m = 1
  for (run in 1:nrun){
    theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
    D_A = matrix(pgdraw(1, theta_A), nrow = n,ncol = n)
    D_A = (D_A + t(D_A)) / 2 
    
    Theta_Y1 = Z %*% B1
    x10 = matrix(-10,n,q1)
    b11_temp = outer(rep(1,n),b11)
    b12_temp = outer(rep(1,n),b12)
    b13_temp = outer(rep(1,n),b13)
    b14_temp = outer(rep(1,n),b14)
    x11 = matrix(10,n,q1)
    phi11<-matrix(rtruncnorm(1, a=x10, b=b11_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#1
    phi11 = ifelse(Y1==1,phi11,0)
    phi12<-matrix(rtruncnorm(1, a=b11_temp, b=b12_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#2
    phi12 = ifelse(Y1==2,phi12,0)
    phi13<-matrix(rtruncnorm(1, a=b12_temp, b=b13_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#3
    phi13 = ifelse(Y1==3,phi13,0)
    phi14<-matrix(rtruncnorm(1, a=b13_temp, b=b14_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#4
    phi14 = ifelse(Y1==4,phi14,0)
    phi15<-matrix(rtruncnorm(1, a=b14_temp, b=x11, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#5
    phi15 = ifelse(Y1==5,phi15,0)
    Phi1 = phi11 + phi12 + phi13 + phi14 + phi15   
    Theta_w1 = Z %*% B1 - Phi1
    W1 = matrix(pgdraw(2,Theta_w1),n,q1) 
    
    Theta_Y2 = Z %*% B2
    x20 = matrix(-10,n,q2)
    b21_temp = outer(rep(1,n),b21)
    b22_temp = outer(rep(1,n),b22)
    x21 = matrix(10,n,q2)
    phi21<-matrix(rtruncnorm(1, a=x20, b=b21_temp, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#1
    phi21 = ifelse(Y2==1,phi21,0)
    phi22<-matrix(rtruncnorm(1, a=b21_temp, b=b22_temp, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#2
    phi22 = ifelse(Y2==2,phi22,0)
    phi23<-matrix(rtruncnorm(1, a=b22_temp, b=x21, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#3
    phi23 = ifelse(Y2==3,phi23,0)
    Phi2 = phi21 + phi22 + phi23
    Theta_w2 = Z %*% B2 - Phi2
    W2 = matrix(pgdraw(2,Theta_w2),n,q2) 
    
    Theta_Y3 = Z %*% B3
    x30 = matrix(-10,n,q3)
    b31_temp = outer(rep(1,n),b31)
    b32_temp = outer(rep(1,n),b32)
    b33_temp = outer(rep(1,n),b33)
    x31 = matrix(10,n,q3)
    phi31<-matrix(rtruncnorm(1, a=x30, b=b31_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#1
    phi31 = ifelse(Y3==1,phi31,0)
    phi32<-matrix(rtruncnorm(1, a=b31_temp, b=b32_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#2
    phi32 = ifelse(Y3==2,phi32,0)
    phi33<-matrix(rtruncnorm(1, a=b32_temp, b=b33_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#3
    phi33 = ifelse(Y3==3,phi33,0)
    phi34<-matrix(rtruncnorm(1, a=b33_temp, b=x31, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#4
    phi34 = ifelse(Y3==4,phi34,0)
    Phi3 = phi31 + phi32 + phi33 + phi34   
    Theta_w3 = Z %*% B3 - Phi3
    W3 = matrix(pgdraw(2,Theta_w3),n,q3)    
    
    theta_Y4 = Z %*% B4 + (matrix(1,n,1) %*% matrix(gamma4,1,q4))
    D_Y4 = matrix(pgdraw(1,theta_Y4),nrow = n,ncol = q4)   
    
    Z_temp = Z %*% t(Z)
    #update alpha
    for(i in 1:n){
      sigma_alphai = 1 / (sum(D_A[i,]) - D_A[i,i] + 1/100)#
      u_temp = A[i,] - 0.5 - D_A[i,] * (alpha + Z_temp[i,])
      u_alphai = sigma_alphai * (sum(u_temp) - u_temp[i])
      alpha[i] = rnorm(1,u_alphai,sqrt(sigma_alphai))
    }
    
    for(i in 1:n){
      Z_i = Z[-i,]
      D_Ai = (D_A[i,])[-i]
      alp_cons = alpha[i] + alpha[-i]
      kappa_Ai = (A[i,])[-i] - 0.5
      
      W1_i = W1[i,]
      phi1_i = Phi1[i,]
      
      W2_i = W2[i,]
      phi2_i = Phi2[i,]
      
      W3_i = W3[i,]
      phi3_i = Phi3[i,]
      
      D_Y4i = D_Y4[i,]
      kappa_Y4i = Y4[i,] - 0.5 
      
      Sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) 
                               + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i
                               + B1 %*% diag(W1_i,nrow = q1) %*% t(B1)
                               + B2 %*% diag(W2_i,nrow = q2) %*% t(B2)
                               + B3 %*% diag(W3_i,nrow = q3) %*% t(B3)
                               + B4 %*% diag(D_Y4i,nrow = q4) %*% t(B4)))
      u_Zi = Sigma_Zi %*% ((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) 
                           + (B1 %*%diag(W1_i,nrow = q1)%*%phi1_i) 
                           + (B2 %*%diag(W2_i,nrow = q2)%*%phi2_i)
                           + (B3 %*%diag(W3_i,nrow = q3)%*%phi3_i)
                           + (B4 %*%(kappa_Y4i - diag(D_Y4i,nrow = q4) %*% gamma4)))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)
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
    
    for (j in 1:q1) {
      W1_j = W1[,j]
      phi1_j = Phi1[,j]
      sigma_B1j = chol2inv(chol(t(Z)%*%diag(W1_j,nrow = n)%*%Z + diag(H, nrow = H)))
      u_B1j = sigma_B1j%*%t(Z)%*%diag(W1_j,nrow = n)%*%phi1_j
      B1[,j] = mvrnorm(1,u_B1j,sigma_B1j)
    }
    
    for (j in 1:q2) {
      W2_j = W2[,j]
      phi2_j = Phi2[,j]
      sigma_B2j = chol2inv(chol(t(Z)%*%diag(W2_j,nrow = n)%*%Z + diag(H, nrow = H)))
      u_B2j = sigma_B2j%*%t(Z)%*%diag(W2_j,nrow = n)%*%phi2_j
      B2[,j] = mvrnorm(1,u_B2j,sigma_B2j)
    }
    
    for (j in 1:q3) {
      W3_j = W3[,j]
      phi3_j = Phi3[,j]
      sigma_B3j = chol2inv(chol(t(Z)%*%diag(W3_j,nrow = n)%*%Z + diag(H, nrow = H)))
      u_B3j = sigma_B3j%*%t(Z)%*%diag(W3_j,nrow = n)%*%phi3_j
      B3[,j] = mvrnorm(1,u_B3j,sigma_B3j)
    }    
    
    for(j in 1:q4){
      D_Y4j = D_Y4[,j]
      sigma_B4j = chol2inv(chol(diag(H, nrow = H) + t(Z) %*% diag(D_Y4j, nrow = n) %*% Z))
      u_B4j = sigma_B4j %*% t(Z) %*% (Y4[,j] - 0.5 - gamma4[j] * D_Y4j)
      B4[,j] = mvrnorm(1, u_B4j, sigma_B4j)
    }
    
    pb1<-ifelse(Y1==1,Phi1,NA)
    pb2<-ifelse(Y1==2,Phi1,NA)
    pb3<-ifelse(Y1==3,Phi1,NA)
    pb4<-ifelse(Y1==4,Phi1,NA)
    pb5<-ifelse(Y1==5,Phi1,NA)
    pb11<-apply(pb1,2,safe_max)
    pb20<-apply(pb2,2,safe_min)
    pb21<-apply(pb2,2,safe_max)
    pb30<-apply(pb3,2,safe_min)
    pb31<-apply(pb3,2,safe_max)
    pb40<-apply(pb4,2,safe_min)
    pb41<-apply(pb4,2,safe_max)
    pb50<-apply(pb5,2,safe_min)
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
    b10 = runif(q1,pb11,pb20)
    b20 = runif(q1,pb21,pb30)
    b30 = runif(q1,pb31,pb40)
    b40 = runif(q1,pb41,pb50)
    b11 = rtruncnorm(1, a=pb11, b=pb20, mean = b10, sd = 10^5)
    b12 = rtruncnorm(1, a=pb21, b=pb30, mean = b20, sd = 10^5)
    b13 = rtruncnorm(1, a=pb31, b=pb40, mean = b30, sd = 10^5)
    b14 = rtruncnorm(1, a=pb41, b=pb50, mean = b40, sd = 10^5)
    
    pb1<-ifelse(Y2==1,Phi2,NA)
    pb2<-ifelse(Y2==2,Phi2,NA)
    pb3<-ifelse(Y2==3,Phi2,NA)
    pb11<-apply(pb1,2,safe_max)
    pb20<-apply(pb2,2,safe_min)
    pb21<-apply(pb2,2,safe_max)
    pb30<-apply(pb3,2,safe_min)
    pb11<-pmax(pb11,-3)
    pb20<-pmin(pb20,0)
    pb21<-pmax(pb21,0)
    pb30<-pmin(pb30,3)
    Q<-matrix(c(pb11,pb20),ncol=2)
    pb11 = apply(Q, 1, min)
    pb20 = apply(Q, 1, max)
    Q<-matrix(c(pb21,pb30),ncol=2)
    pb21 = apply(Q, 1, min)
    pb30 = apply(Q, 1, max)
    b10 = runif(q2,pb11,pb20)
    b20 = runif(q2,pb21,pb30)
    b21 = rtruncnorm(1, a=pb11, b=pb20, mean = b10, sd = 10^5)
    b22 = rtruncnorm(1, a=pb21, b=pb30, mean = b20, sd = 10^5)
    
    pb1<-ifelse(Y3==1,Phi3,NA)
    pb2<-ifelse(Y3==2,Phi3,NA)
    pb3<-ifelse(Y3==3,Phi3,NA)
    pb4<-ifelse(Y3==4,Phi3,NA)
    pb11<-apply(pb1,2,safe_max)
    pb20<-apply(pb2,2,safe_min)
    pb21<-apply(pb2,2,safe_max)
    pb30<-apply(pb3,2,safe_min)
    pb31<-apply(pb3,2,safe_max)
    pb40<-apply(pb4,2,safe_min)
    pb11<-pmax(pb11,-3)
    pb20<-pmin(pb20,-1)
    pb21<-pmax(pb21,-1)
    pb30<-pmin(pb30,1)
    pb31<-pmax(pb31,1)
    pb40<-pmin(pb40,3)
    Q<-matrix(c(pb11,pb20),ncol=2)
    pb11 = apply(Q, 1, min)
    pb20 = apply(Q, 1, max)
    Q<-matrix(c(pb21,pb30),ncol=2)
    pb21 = apply(Q, 1, min)
    pb30 = apply(Q, 1, max)
    Q<-matrix(c(pb31,pb40),ncol=2)
    pb31 = apply(Q, 1, min)
    pb40 = apply(Q, 1, max)
    b10 = runif(q3,pb11,pb20)
    b20 = runif(q3,pb21,pb30)
    b30 = runif(q3,pb31,pb40)
    b31 = rtruncnorm(1, a=pb11, b=pb20, mean = b10, sd = 10^5)
    b32 = rtruncnorm(1, a=pb21, b=pb30, mean = b20, sd = 10^5)
    b33 = rtruncnorm(1, a=pb31, b=pb40, mean = b30, sd = 10^5)
    
    si_g4 = apply(D_Y4, 2, sum) + 1/100
    sigma_gamma4 = diag(1 / si_g4)
    u_gamma4_temp = Y4 - 0.5 - D_Y4 * (Z %*% B4)
    u_gamma4 = sigma_gamma4 %*% (apply(u_gamma4_temp, 2, sum))
    gamma4 = mvrnorm(1, u_gamma4, sigma_gamma4)
    
    Y1_hat = get_item_Y1(Z,B1,b11,b12,b13,b14)
    Y1[MissY1 == 1] = Y1_hat[MissY1 == 1]
    
    Y3_hat = get_item_Y3(Z,B3,b31,b32,b33)
    Y3[MissY3 == 1] = Y3_hat[MissY3 == 1]
    
    active = which(zta > c(1:H))
    Hstar = length(active)
    if (run >= start_adapt & u[run] <= exp(alpha0 + alpha1 * run)){
      if (Hstar < H - 1){
        # set truncation to Hstar[t] and subset all variables, keeping only active columns
        H = Hstar + 1
        theta_inv = c(theta_inv[active],theta_inf^(-1))
        w = c(w[active],1-sum(w[active]))
        Z = cbind(Z[, active, drop = FALSE], matrix(rnorm(n, mean=0, sd=sqrt(theta_inf)), ncol = 1))
        B1 = rbind(B1[active, , drop = FALSE], matrix(0.05 * rnorm(q1), nrow = 1))
        B2 = rbind(B2[active, , drop = FALSE], matrix(0.05 * rnorm(q2), nrow = 1))
        B3 = rbind(B3[active, , drop = FALSE], matrix(0.05 * rnorm(q3), nrow = 1))
        B4 = rbind(B4[active, , drop = FALSE], matrix(0.05 * rnorm(q4), nrow = 1))
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
        Z = cbind(Z, matrix(rnorm(n, mean=0, sd=sqrt(theta_inf)), ncol = 1))
        B1 = rbind(B1, matrix(0.05 * rnorm(q1), nrow = 1))
        B2 = rbind(B2, matrix(0.05 * rnorm(q2), nrow = 1))
        B3 = rbind(B3, matrix(0.05 * rnorm(q3), nrow = 1))
        B4 = rbind(B4, matrix(0.05 * rnorm(q4), nrow = 1))
        #theta_inv_B = c(theta_inv_B,rgamma(1,a_theta_B,b_theta_B))
        #B = rbind(B,rnorm(q,0,sqrt(theta_inv_B[H])))
        zta = c(zta,H-1)
      }
    }
    
    H_hat[run] = Hstar
    if((run > burn) &((run-burn) %% thin == 0)){
      b11_hat[m,] = b11
      b12_hat[m,] = b12
      b13_hat[m,] = b13
      b14_hat[m,] = b14
      b21_hat[m,] = b21
      b22_hat[m,] = b22
      b31_hat[m,] = b31
      b32_hat[m,] = b32
      b33_hat[m,] = b33
      gamma4_hat[m,] = gamma4
      if(Hstar>0){
        B1_hat[[m]] = B1[1:Hstar,, drop=FALSE]
        B2_hat[[m]] = B2[1:Hstar,, drop=FALSE]
        B3_hat[[m]] = B3[1:Hstar,, drop=FALSE]
        B4_hat[[m]] = B4[1:Hstar,, drop=FALSE]
        Z_hat[[m]] = Z[,1:Hstar, drop=FALSE]
      }else{
        B1_hat[[m]] = B1
        B2_hat[[m]] = B2
        B3_hat[[m]] = B3
        B4_hat[[m]] = B4
        Z_hat[[m]] = Z
      }
      alpha_hat[m,] = alpha
      m = m+1
    }    
  }
  output = list("H" = H_hat,"alpha" = alpha_hat,"B1" = B1_hat,"B2" = B2_hat,"B3" = B3_hat,"B4" = B4_hat,
                "b11" = b11_hat,"b12" = b12_hat,"b13" = b13_hat,"b14" = b14_hat,
                "b21" = b21_hat,"b22" = b22_hat,"b31" = b31_hat,"b32" = b32_hat,"b33" = b33_hat,
                "gamma4" = gamma4_hat,"Z" = Z_hat)
  return(output)
}


network_mix_real_data_kf = function(A,Y1,Y2,Y3,Y4,MissY1,MissY3,H,nrun,burn,thin){
  n = dim(Y1)[1]
  q1 = dim(Y1)[2]
  q2 = dim(Y2)[2]
  q3 = dim(Y3)[2]
  q4 = dim(Y4)[2]
  
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)  
  theta_inv = rep(1,H)  
  
  theta_inv_B1 = rep(1,H)
  theta_inv_B2 = rep(1,H)
  theta_inv_B3 = rep(1,H)
  theta_inv_B4 = rep(1,H)
  w = rep(1,H)
  zta = rep(1,H)    
  B1 = matrix(rnorm(H * q1),nrow = H, ncol = q1)
  B2 = matrix(rnorm(H * q2),nrow = H, ncol = q2)
  B3 = matrix(rnorm(H * q3),nrow = H, ncol = q3)
  B4 = matrix(rnorm(H * q4),nrow = H, ncol = q4)
  
  Phi1 = matrix(1,n,q1)#
  W1 = matrix(1,n,q1)
  b11<-rtruncnorm(q1, a=-3, b=-1.5, mean = 0, sd = 1)
  b12<-rtruncnorm(q1, a=-1.5, b=0, mean = 0, sd = 1)
  b13<-rtruncnorm(q1, a=0, b=1.5, mean = 0, sd = 1)
  b14<-rtruncnorm(q1, a=1.5, b=3, mean = 0, sd = 1)  
  
  Phi2 = matrix(1,n,q2)#
  W2 = matrix(1,n,q2)
  b21<-rtruncnorm(q2, a=-3, b=0, mean = 0, sd = 1)
  b22<-rtruncnorm(q2, a=0, b=3, mean = 0, sd = 1)
  
  Phi3 = matrix(1,n,q3)#
  W3 = matrix(1,n,q3)
  b31<-rtruncnorm(q3, a=-3, b=-1, mean = 0, sd = 1)
  b32<-rtruncnorm(q3, a=-1, b=1, mean = 0, sd = 1)
  b33<-rtruncnorm(q3, a=1, b=3, mean = 0, sd = 1)
  
  gamma4 = rep(1,q4)
  
  N_sample = ceiling((nrun - burn)/thin)
  
  Z_hat = array(0,dim = c(N_sample,n,H))
  B1_hat = array(0,dim = c(N_sample,H,q1))
  B2_hat = array(0,dim = c(N_sample,H,q2))
  B3_hat = array(0,dim = c(N_sample,H,q3))
  B4_hat = array(0,dim = c(N_sample,H,q4))
  
  alpha_hat = matrix(0,N_sample,n)  
  
  b11_hat = matrix(0,N_sample,q1)
  b12_hat = matrix(0,N_sample,q1)
  b13_hat = matrix(0,N_sample,q1)
  b14_hat = matrix(0,N_sample,q1)
  
  b21_hat = matrix(0,N_sample,q2)
  b22_hat = matrix(0,N_sample,q2)
  
  b31_hat = matrix(0,N_sample,q3)
  b32_hat = matrix(0,N_sample,q3)
  b33_hat = matrix(0,N_sample,q3)
  
  gamma4_hat = matrix(0,N_sample,q4)  
  
  m = 1
  for (run in 1:nrun){
    theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
    D_A = matrix(pgdraw(1, theta_A), nrow = n,ncol = n)
    D_A = (D_A + t(D_A)) / 2 
    
    Theta_Y1 = Z %*% B1
    x10 = matrix(-10,n,q1)
    b11_temp = outer(rep(1,n),b11)
    b12_temp = outer(rep(1,n),b12)
    b13_temp = outer(rep(1,n),b13)
    b14_temp = outer(rep(1,n),b14)
    x11 = matrix(10,n,q1)
    phi11<-matrix(rtruncnorm(1, a=x10, b=b11_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#1
    phi11 = ifelse(Y1==1,phi11,0)
    phi12<-matrix(rtruncnorm(1, a=b11_temp, b=b12_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#2
    phi12 = ifelse(Y1==2,phi12,0)
    phi13<-matrix(rtruncnorm(1, a=b12_temp, b=b13_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#3
    phi13 = ifelse(Y1==3,phi13,0)
    phi14<-matrix(rtruncnorm(1, a=b13_temp, b=b14_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#4
    phi14 = ifelse(Y1==4,phi14,0)
    phi15<-matrix(rtruncnorm(1, a=b14_temp, b=x11, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#5
    phi15 = ifelse(Y1==5,phi15,0)
    Phi1 = phi11 + phi12 + phi13 + phi14 + phi15   
    Theta_w1 = Z %*% B1 - Phi1
    W1 = matrix(pgdraw(2,Theta_w1),n,q1) 
    
    Theta_Y2 = Z %*% B2
    x20 = matrix(-10,n,q2)
    b21_temp = outer(rep(1,n),b21)
    b22_temp = outer(rep(1,n),b22)
    x21 = matrix(10,n,q2)
    phi21<-matrix(rtruncnorm(1, a=x20, b=b21_temp, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#1
    phi21 = ifelse(Y2==1,phi21,0)
    phi22<-matrix(rtruncnorm(1, a=b21_temp, b=b22_temp, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#2
    phi22 = ifelse(Y2==2,phi22,0)
    phi23<-matrix(rtruncnorm(1, a=b22_temp, b=x21, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#3
    phi23 = ifelse(Y2==3,phi23,0)
    Phi2 = phi21 + phi22 + phi23
    Theta_w2 = Z %*% B2 - Phi2
    W2 = matrix(pgdraw(2,Theta_w2),n,q2) 
    
    Theta_Y3 = Z %*% B3
    x30 = matrix(-10,n,q3)
    b31_temp = outer(rep(1,n),b31)
    b32_temp = outer(rep(1,n),b32)
    b33_temp = outer(rep(1,n),b33)
    x31 = matrix(10,n,q3)
    phi31<-matrix(rtruncnorm(1, a=x30, b=b31_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#1
    phi31 = ifelse(Y3==1,phi31,0)
    phi32<-matrix(rtruncnorm(1, a=b31_temp, b=b32_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#2
    phi32 = ifelse(Y3==2,phi32,0)
    phi33<-matrix(rtruncnorm(1, a=b32_temp, b=b33_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#3
    phi33 = ifelse(Y3==3,phi33,0)
    phi34<-matrix(rtruncnorm(1, a=b33_temp, b=x31, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#4
    phi34 = ifelse(Y3==4,phi34,0)
    Phi3 = phi31 + phi32 + phi33 + phi34   
    Theta_w3 = Z %*% B3 - Phi3
    W3 = matrix(pgdraw(2,Theta_w3),n,q3)    
    
    theta_Y4 = Z %*% B4 + (matrix(1,n,1) %*% matrix(gamma4,1,q4))
    D_Y4 = matrix(pgdraw(1,theta_Y4),nrow = n,ncol = q4)   
    
    Z_temp = Z %*% t(Z)
    #update alpha
    for(i in 1:n){
      sigma_alphai = 1 / (sum(D_A[i,]) - D_A[i,i] + 1/100)#
      u_temp = A[i,] - 0.5 - D_A[i,] * (alpha + Z_temp[i,])
      u_alphai = sigma_alphai * (sum(u_temp) - u_temp[i])
      alpha[i] = rnorm(1,u_alphai,sqrt(sigma_alphai))
    }
    
    for(i in 1:n){
      Z_i = Z[-i,]
      D_Ai = (D_A[i,])[-i]
      alp_cons = alpha[i] + alpha[-i]
      kappa_Ai = (A[i,])[-i] - 0.5
      
      W1_i = W1[i,]
      phi1_i = Phi1[i,]
      
      W2_i = W2[i,]
      phi2_i = Phi2[i,]
      
      W3_i = W3[i,]
      phi3_i = Phi3[i,]
      
      D_Y4i = D_Y4[i,]
      kappa_Y4i = Y4[i,] - 0.5 
      
      Sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) 
                               + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i
                               + B1 %*% diag(W1_i,nrow = q1) %*% t(B1)
                               + B2 %*% diag(W2_i,nrow = q2) %*% t(B2)
                               + B3 %*% diag(W3_i,nrow = q3) %*% t(B3)
                               + B4 %*% diag(D_Y4i,nrow = q4) %*% t(B4)))
      u_Zi = Sigma_Zi %*% ((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) 
                           + (B1 %*%diag(W1_i,nrow = q1)%*%phi1_i) 
                           + (B2 %*%diag(W2_i,nrow = q2)%*%phi2_i)
                           + (B3 %*%diag(W3_i,nrow = q3)%*%phi3_i)
                           + (B4 %*%(kappa_Y4i - diag(D_Y4i,nrow = q4) %*% gamma4)))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)
    }
    
    for (j in 1:q1) {
      W1_j = W1[,j]
      phi1_j = Phi1[,j]
      sigma_B1j = chol2inv(chol(t(Z)%*%diag(W1_j,nrow = n)%*%Z + diag(H, nrow = H)))
      u_B1j = sigma_B1j%*%t(Z)%*%diag(W1_j,nrow = n)%*%phi1_j
      B1[,j] = mvrnorm(1,u_B1j,sigma_B1j)
    }
    
    for (j in 1:q2) {
      W2_j = W2[,j]
      phi2_j = Phi2[,j]
      sigma_B2j = chol2inv(chol(t(Z)%*%diag(W2_j,nrow = n)%*%Z + diag(H, nrow = H)))
      u_B2j = sigma_B2j%*%t(Z)%*%diag(W2_j,nrow = n)%*%phi2_j
      B2[,j] = mvrnorm(1,u_B2j,sigma_B2j)
    }
    
    for (j in 1:q3) {
      W3_j = W3[,j]
      phi3_j = Phi3[,j]
      sigma_B3j = chol2inv(chol(t(Z)%*%diag(W3_j,nrow = n)%*%Z + diag(H, nrow = H)))
      u_B3j = sigma_B3j%*%t(Z)%*%diag(W3_j,nrow = n)%*%phi3_j
      B3[,j] = mvrnorm(1,u_B3j,sigma_B3j)
    }    
    
    for(j in 1:q4){
      D_Y4j = D_Y4[,j]
      sigma_B4j = chol2inv(chol(diag(H, nrow = H) + t(Z) %*% diag(D_Y4j, nrow = n) %*% Z))
      u_B4j = sigma_B4j %*% t(Z) %*% (Y4[,j] - 0.5 - gamma4[j] * D_Y4j)
      B4[,j] = mvrnorm(1, u_B4j, sigma_B4j)
    }
    
    pb1<-ifelse(Y1==1,Phi1,NA)
    pb2<-ifelse(Y1==2,Phi1,NA)
    pb3<-ifelse(Y1==3,Phi1,NA)
    pb4<-ifelse(Y1==4,Phi1,NA)
    pb5<-ifelse(Y1==5,Phi1,NA)
    pb11<-apply(pb1,2,safe_max)
    pb20<-apply(pb2,2,safe_min)
    pb21<-apply(pb2,2,safe_max)
    pb30<-apply(pb3,2,safe_min)
    pb31<-apply(pb3,2,safe_max)
    pb40<-apply(pb4,2,safe_min)
    pb41<-apply(pb4,2,safe_max)
    pb50<-apply(pb5,2,safe_min)
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
    b10 = runif(q1,pb11,pb20)
    b20 = runif(q1,pb21,pb30)
    b30 = runif(q1,pb31,pb40)
    b40 = runif(q1,pb41,pb50)
    b11 = rtruncnorm(1, a=pb11, b=pb20, mean = b10, sd = 10^5)
    b12 = rtruncnorm(1, a=pb21, b=pb30, mean = b20, sd = 10^5)
    b13 = rtruncnorm(1, a=pb31, b=pb40, mean = b30, sd = 10^5)
    b14 = rtruncnorm(1, a=pb41, b=pb50, mean = b40, sd = 10^5)
    
    pb1<-ifelse(Y2==1,Phi2,NA)
    pb2<-ifelse(Y2==2,Phi2,NA)
    pb3<-ifelse(Y2==3,Phi2,NA)
    pb11<-apply(pb1,2,safe_max)
    pb20<-apply(pb2,2,safe_min)
    pb21<-apply(pb2,2,safe_max)
    pb30<-apply(pb3,2,safe_min)
    pb11<-pmax(pb11,-3)
    pb20<-pmin(pb20,0)
    pb21<-pmax(pb21,0)
    pb30<-pmin(pb30,3)
    Q<-matrix(c(pb11,pb20),ncol=2)
    pb11 = apply(Q, 1, min)
    pb20 = apply(Q, 1, max)
    Q<-matrix(c(pb21,pb30),ncol=2)
    pb21 = apply(Q, 1, min)
    pb30 = apply(Q, 1, max)
    b10 = runif(q2,pb11,pb20)
    b20 = runif(q2,pb21,pb30)
    b21 = rtruncnorm(1, a=pb11, b=pb20, mean = b10, sd = 10^5)
    b22 = rtruncnorm(1, a=pb21, b=pb30, mean = b20, sd = 10^5)
    
    pb1<-ifelse(Y3==1,Phi3,NA)
    pb2<-ifelse(Y3==2,Phi3,NA)
    pb3<-ifelse(Y3==3,Phi3,NA)
    pb4<-ifelse(Y3==4,Phi3,NA)
    pb11<-apply(pb1,2,safe_max)
    pb20<-apply(pb2,2,safe_min)
    pb21<-apply(pb2,2,safe_max)
    pb30<-apply(pb3,2,safe_min)
    pb31<-apply(pb3,2,safe_max)
    pb40<-apply(pb4,2,safe_min)
    pb11<-pmax(pb11,-3)
    pb20<-pmin(pb20,-1)
    pb21<-pmax(pb21,-1)
    pb30<-pmin(pb30,1)
    pb31<-pmax(pb31,1)
    pb40<-pmin(pb40,3)
    Q<-matrix(c(pb11,pb20),ncol=2)
    pb11 = apply(Q, 1, min)
    pb20 = apply(Q, 1, max)
    Q<-matrix(c(pb21,pb30),ncol=2)
    pb21 = apply(Q, 1, min)
    pb30 = apply(Q, 1, max)
    Q<-matrix(c(pb31,pb40),ncol=2)
    pb31 = apply(Q, 1, min)
    pb40 = apply(Q, 1, max)
    b10 = runif(q3,pb11,pb20)
    b20 = runif(q3,pb21,pb30)
    b30 = runif(q3,pb31,pb40)
    b31 = rtruncnorm(1, a=pb11, b=pb20, mean = b10, sd = 10^5)
    b32 = rtruncnorm(1, a=pb21, b=pb30, mean = b20, sd = 10^5)
    b33 = rtruncnorm(1, a=pb31, b=pb40, mean = b30, sd = 10^5)
    
    si_g4 = apply(D_Y4, 2, sum) + 1/100
    sigma_gamma4 = diag(1 / si_g4)
    u_gamma4_temp = Y4 - 0.5 - D_Y4 * (Z %*% B4)
    u_gamma4 = sigma_gamma4 %*% (apply(u_gamma4_temp, 2, sum))
    gamma4 = mvrnorm(1, u_gamma4, sigma_gamma4)
    
    Y1_hat = get_item_Y1(Z,B1,b11,b12,b13,b14)
    Y1[MissY1 == 1] = Y1_hat[MissY1 == 1]
    
    Y3_hat = get_item_Y3(Z,B3,b31,b32,b33)
    Y3[MissY3 == 1] = Y3_hat[MissY3 == 1]  
    
    if((run > burn) &((run-burn) %% thin == 0)){
      b11_hat[m,] = b11
      b12_hat[m,] = b12
      b13_hat[m,] = b13
      b14_hat[m,] = b14
      b21_hat[m,] = b21
      b22_hat[m,] = b22
      b31_hat[m,] = b31
      b32_hat[m,] = b32
      b33_hat[m,] = b33
      gamma4_hat[m,] = gamma4
      
      Z_hat[m,,] = Z
      B1_hat[m,,] = B1    
      B2_hat[m,,] = B2  
      B3_hat[m,,] = B3  
      B4_hat[m,,] = B4  
      
      alpha_hat[m,] = alpha
      m = m+1      
    }    
  }
  output = list("alpha" = alpha_hat,"B1" = B1_hat,"B2" = B2_hat,"B3" = B3_hat,"B4" = B4_hat,
                "b11" = b11_hat,"b12" = b12_hat,"b13" = b13_hat,"b14" = b14_hat,
                "b21" = b21_hat,"b22" = b22_hat,"b31" = b31_hat,"b32" = b32_hat,"b33" = b33_hat,
                "gamma4" = gamma4_hat,"Z" = Z_hat)
  return(output)
}

network_mix_real_data_kf_getZ = function(A,Y1,Y2,Y3,Y4,MissY1,MissY3,B1,B2,B3,B4,b11,b12,b13,
                                         b14,b21,b22,b31,b32,b33,gamma4,H,nrun,burn,thin){
  n = dim(Y1)[1]
  q1 = dim(Y1)[2]
  q2 = dim(Y2)[2]
  q3 = dim(Y3)[2]
  q4 = dim(Y4)[2]
  
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)  
  theta_inv = rep(1,H)  
  
  theta_inv_B1 = rep(1,H)
  theta_inv_B2 = rep(1,H)
  theta_inv_B3 = rep(1,H)
  theta_inv_B4 = rep(1,H)
  w = rep(1,H)
  zta = rep(1,H)  
  
  Phi1 = matrix(1,n,q1)#文章中的Z
  W1 = matrix(1,n,q1)
  
  Phi2 = matrix(1,n,q2)#文章中的Z
  W2 = matrix(1,n,q2)
  
  Phi3 = matrix(1,n,q3)#文章中的Z
  W3 = matrix(1,n,q3)
  N_sample = ceiling((nrun - burn)/thin)
  
  Z_hat = array(0,dim = c(N_sample,n,H))
  alpha_hat = matrix(0,N_sample,n)    
  m = 1
  for (run in 1:nrun){
    theta_A = Z %*% t(Z) + matrix(1,n,1)%*%matrix(alpha,1,n) + matrix(alpha,n,1) %*%matrix(1,1,n)
    D_A = matrix(pgdraw(1, theta_A), nrow = n,ncol = n)
    D_A = (D_A + t(D_A)) / 2 
    
    Theta_Y1 = Z %*% B1
    x10 = matrix(-10,n,q1)
    b11_temp = outer(rep(1,n),b11)
    b12_temp = outer(rep(1,n),b12)
    b13_temp = outer(rep(1,n),b13)
    b14_temp = outer(rep(1,n),b14)
    x11 = matrix(10,n,q1)
    phi11<-matrix(rtruncnorm(1, a=x10, b=b11_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#1
    phi11 = ifelse(Y1==1,phi11,0)
    phi12<-matrix(rtruncnorm(1, a=b11_temp, b=b12_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#2
    phi12 = ifelse(Y1==2,phi12,0)
    phi13<-matrix(rtruncnorm(1, a=b12_temp, b=b13_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#3
    phi13 = ifelse(Y1==3,phi13,0)
    phi14<-matrix(rtruncnorm(1, a=b13_temp, b=b14_temp, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#4
    phi14 = ifelse(Y1==4,phi14,0)
    phi15<-matrix(rtruncnorm(1, a=b14_temp, b=x11, mean = Theta_Y1, sd = sqrt(1/W1)),n,q1)#5
    phi15 = ifelse(Y1==5,phi15,0)
    Phi1 = phi11 + phi12 + phi13 + phi14 + phi15   
    Theta_w1 = Z %*% B1 - Phi1
    W1 = matrix(pgdraw(2,Theta_w1),n,q1) 
    
    Theta_Y2 = Z %*% B2
    x20 = matrix(-10,n,q2)
    b21_temp = outer(rep(1,n),b21)
    b22_temp = outer(rep(1,n),b22)
    x21 = matrix(10,n,q2)
    phi21<-matrix(rtruncnorm(1, a=x20, b=b21_temp, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#1
    phi21 = ifelse(Y2==1,phi21,0)
    phi22<-matrix(rtruncnorm(1, a=b21_temp, b=b22_temp, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#2
    phi22 = ifelse(Y2==2,phi22,0)
    phi23<-matrix(rtruncnorm(1, a=b22_temp, b=x21, mean = Theta_Y2, sd = sqrt(1/W2)),n,q2)#3
    phi23 = ifelse(Y2==3,phi23,0)
    Phi2 = phi21 + phi22 + phi23
    Theta_w2 = Z %*% B2 - Phi2
    W2 = matrix(pgdraw(2,Theta_w2),n,q2) 
    
    Theta_Y3 = Z %*% B3
    x30 = matrix(-10,n,q3)
    b31_temp = outer(rep(1,n),b31)
    b32_temp = outer(rep(1,n),b32)
    b33_temp = outer(rep(1,n),b33)
    x31 = matrix(10,n,q3)
    phi31<-matrix(rtruncnorm(1, a=x30, b=b31_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#1
    phi31 = ifelse(Y3==1,phi31,0)
    phi32<-matrix(rtruncnorm(1, a=b31_temp, b=b32_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#2
    phi32 = ifelse(Y3==2,phi32,0)
    phi33<-matrix(rtruncnorm(1, a=b32_temp, b=b33_temp, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#3
    phi33 = ifelse(Y3==3,phi33,0)
    phi34<-matrix(rtruncnorm(1, a=b33_temp, b=x31, mean = Theta_Y3, sd = sqrt(1/W3)),n,q3)#4
    phi34 = ifelse(Y3==4,phi34,0)
    Phi3 = phi31 + phi32 + phi33 + phi34   
    Theta_w3 = Z %*% B3 - Phi3
    W3 = matrix(pgdraw(2,Theta_w3),n,q3)    
    
    theta_Y4 = Z %*% B4 + (matrix(1,n,1) %*% matrix(gamma4,1,q4))
    D_Y4 = matrix(pgdraw(1,theta_Y4),nrow = n,ncol = q4)   
    
    Z_temp = Z %*% t(Z)
    #update alpha
    for(i in 1:n){
      sigma_alphai = 1 / (sum(D_A[i,]) - D_A[i,i] + 1/100)#
      u_temp = A[i,] - 0.5 - D_A[i,] * (alpha + Z_temp[i,])
      u_alphai = sigma_alphai * (sum(u_temp) - u_temp[i])
      alpha[i] = rnorm(1,u_alphai,sqrt(sigma_alphai))
    }
    
    for(i in 1:n){
      Z_i = Z[-i,]
      D_Ai = (D_A[i,])[-i]
      alp_cons = alpha[i] + alpha[-i]
      kappa_Ai = (A[i,])[-i] - 0.5
      
      W1_i = W1[i,]
      phi1_i = Phi1[i,]
      
      W2_i = W2[i,]
      phi2_i = Phi2[i,]
      
      W3_i = W3[i,]
      phi3_i = Phi3[i,]
      
      D_Y4i = D_Y4[i,]
      kappa_Y4i = Y4[i,] - 0.5 
      
      Sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) 
                               + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i
                               + B1 %*% diag(W1_i,nrow = q1) %*% t(B1)
                               + B2 %*% diag(W2_i,nrow = q2) %*% t(B2)
                               + B3 %*% diag(W3_i,nrow = q3) %*% t(B3)
                               + B4 %*% diag(D_Y4i,nrow = q4) %*% t(B4)))
      u_Zi = Sigma_Zi %*% ((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) 
                           + (B1 %*%diag(W1_i,nrow = q1)%*%phi1_i) 
                           + (B2 %*%diag(W2_i,nrow = q2)%*%phi2_i)
                           + (B3 %*%diag(W3_i,nrow = q3)%*%phi3_i)
                           + (B4 %*%(kappa_Y4i - diag(D_Y4i,nrow = q4) %*% gamma4)))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)
    }
    Y1_hat = get_item_Y1(Z,B1,b11,b12,b13,b14)
    Y1[MissY1 == 1] = Y1_hat[MissY1 == 1]
    
    Y3_hat = get_item_Y3(Z,B3,b31,b32,b33)
    Y3[MissY3 == 1] = Y3_hat[MissY3 == 1]
    
    if((run > burn) &((run-burn) %% thin == 0)){
      Z_hat[m,,] = Z
      alpha_hat[m,] = alpha
      m = m+1      
    } 
  }
  output = list("alpha" = alpha_hat,"Z" = Z_hat)
  return(output)  
}



get_loglikehood_A = function(A,alpha,Z){
  n = nrow(A)
  theta_A = Z %*% t(Z) + outer(alpha, alpha, "+")
  L_A = A * theta_A - log(1 + exp(theta_A))
  ll_A = sum(L_A[upper.tri(L_A)])
  return(ll_A)
}

get_loglikehood_Y1 = function(Y,MissY,b1,b2,b3,b4,B,Z){
  n = dim(Y)[1]
  q = dim(Y)[2]
  XB = Z %*% B 
  cp1 = plogis(XB + matrix(b1, n, q, byrow = TRUE))
  cp2 = plogis(XB + matrix(b2, n, q, byrow = TRUE))
  cp3 = plogis(XB + matrix(b3, n, q, byrow = TRUE))
  cp4 = plogis(XB + matrix(b4, n, q, byrow = TRUE))
  PY = list()
  PY[[1]] = cp1
  PY[[2]] = cp2 - cp1
  PY[[3]] = cp3 - cp2
  PY[[4]] = cp4 - cp3
  PY[[5]] = 1 - cp4
  L_Y_mat = matrix(0, n, q)
  for (k in 1:5) {
    idx = (Y == k)
    L_Y_mat[idx] = PY[[k]][idx]
  }
  L_Y_mat[L_Y_mat < 1e-15] = 1e-15
  L_Y_mat[MissY == 1] = 1
  return(sum(log(L_Y_mat)))
}

get_loglikehood_Y2 = function(Y,b1,b2,B,Z){
  n = dim(Y)[1]
  q = dim(Y)[2]
  XB = Z %*% B 
  cp1 = plogis(XB + matrix(b1, n, q, byrow = TRUE))
  cp2 = plogis(XB + matrix(b2, n, q, byrow = TRUE))
  PY = list()
  PY[[1]] = cp1
  PY[[2]] = cp2 - cp1
  PY[[3]] = 1 - cp2
  L_Y_mat = matrix(0, n, q)
  for (k in 1:3) {
    idx = (Y == k)
    L_Y_mat[idx] = PY[[k]][idx]
  }
  L_Y_mat[L_Y_mat < 1e-15] = 1e-15
  return(sum(log(L_Y_mat)))
}

get_loglikehood_Y3 = function(Y,MissY,b1,b2,b3,B,Z){
  n = dim(Y)[1]
  q = dim(Y)[2]
  XB = Z %*% B 
  cp1 = plogis(XB + matrix(b1, n, q, byrow = TRUE))
  cp2 = plogis(XB + matrix(b2, n, q, byrow = TRUE))
  cp3 = plogis(XB + matrix(b3, n, q, byrow = TRUE))
  PY = list()
  PY[[1]] = cp1
  PY[[2]] = cp2 - cp1
  PY[[3]] = cp3 - cp2
  PY[[4]] = 1 - cp3
  L_Y_mat = matrix(0, n, q)
  for (k in 1:4) {
    idx = (Y == k)
    L_Y_mat[idx] = PY[[k]][idx]
  }
  L_Y_mat[L_Y_mat < 1e-15] = 1e-15
  L_Y_mat[MissY == 1] = 1
  return(sum(log(L_Y_mat)))
}

get_loglikehood_Y4 = function(Y,gamma_Y,B,Z){
  n = dim(Y)[1]
  q = dim(Y)[2]  
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_Y = Y * theta_Y - log(1 + exp(theta_Y))
  return(sum(L_Y))
}

get_loglikehood_data3 = function(A,Y1,Y2,Y3,Y4,MissY1,MissY3,b11,b12,b13,b14,b21,
                                 b22,b31,b32,b33,gamma4,B1,B2,B3,B4,alpha,Z){
  L_A = get_loglikehood_A(A,alpha,Z)
  L_Y1 = get_loglikehood_Y1(Y1,MissY1,b11,b12,b13,b14,B1,Z)
  L_Y2 = get_loglikehood_Y2(Y2,b21,b22,B2,Z)
  L_Y3 = get_loglikehood_Y3(Y3,MissY3,b31,b32,b33,B3,Z)
  L_Y4 = get_loglikehood_Y4(Y4,gamma4,B4,Z)
  return(L_A + L_Y1 + L_Y2 + L_Y3 + L_Y4)
}

get_AIC_data3 = function(A,Y1,Y2,Y3,Y4,MissY1,MissY3,b11,b12,b13,b14,b21,
                         b22,b31,b32,b33,gamma4,B1,B2,B3,B4,alpha,Z){
  n = dim(Y1)[1]
  q1 = dim(Y1)[2]
  q2 = dim(Y2)[2]
  q3 = dim(Y3)[2]
  q4 = dim(Y4)[2]
  k = dim(Z)[2]
  L_AY = get_loglikehood_data3(A,Y1,Y2,Y3,Y4,MissY1,MissY3,b11,b12,b13,b14,b21,
                               b22,b31,b32,b33,gamma4,B1,B2,B3,B4,alpha,Z)
  return(2*(n*(k+1)+(q1+q2+q3+q4)*k+4*q1+2*q2+3*q3+q4)-2*L_AY)
}

get_BIC_data3 = function(A,Y1,Y2,Y3,Y4,MissY1,MissY3,b11,b12,b13,b14,b21,
                         b22,b31,b32,b33,gamma4,B1,B2,B3,B4,alpha,Z){
  n = dim(Y1)[1]
  q1 = dim(Y1)[2]
  q2 = dim(Y2)[2]
  q3 = dim(Y3)[2]
  q4 = dim(Y4)[2]
  k = dim(Z)[2]
  L_AY = get_loglikehood_data3(A,Y1,Y2,Y3,Y4,MissY1,MissY3,b11,b12,b13,b14,b21,
                               b22,b31,b32,b33,gamma4,B1,B2,B3,B4,alpha,Z)
  return(2*(n*(k+1)+(q1+q2+q3+q4)*k+4*q1+2*q2+3*q3+q4)*log(n)-2*L_AY)
}


get_loglikehood_point_A = function(A,alpha,Z){
  n = nrow(A)
  theta_A = Z %*% t(Z) + outer(alpha, alpha, "+")
  L_A = A * theta_A - log(1 + exp(theta_A))
  return(L_A)
}

get_loglikehood_point_Y1 = function(Y,MissY,b1,b2,b3,b4,B,Z){
  n = dim(Y)[1]
  q = dim(Y)[2]
  XB = Z %*% B 
  cp1 = plogis(XB + matrix(b1, n, q, byrow = TRUE))
  cp2 = plogis(XB + matrix(b2, n, q, byrow = TRUE))
  cp3 = plogis(XB + matrix(b3, n, q, byrow = TRUE))
  cp4 = plogis(XB + matrix(b4, n, q, byrow = TRUE))
  PY = list()
  PY[[1]] = cp1
  PY[[2]] = cp2 - cp1
  PY[[3]] = cp3 - cp2
  PY[[4]] = cp4 - cp3
  PY[[5]] = 1 - cp4
  L_Y_mat = matrix(0, n, q)
  for (k in 1:5) {
    idx = (Y == k)
    L_Y_mat[idx] = PY[[k]][idx]
  }
  L_Y_mat[L_Y_mat < 1e-15] = 1e-15
  L_Y_mat[MissY == 1] = 0
  return((log(L_Y_mat)))
}

get_loglikehood_point_Y2 = function(Y,b1,b2,B,Z){
  n = dim(Y)[1]
  q = dim(Y)[2]
  XB = Z %*% B 
  cp1 = plogis(XB + matrix(b1, n, q, byrow = TRUE))
  cp2 = plogis(XB + matrix(b2, n, q, byrow = TRUE))
  PY = list()
  PY[[1]] = cp1
  PY[[2]] = cp2 - cp1
  PY[[3]] = 1 - cp2
  L_Y_mat = matrix(0, n, q)
  for (k in 1:3) {
    idx = (Y == k)
    L_Y_mat[idx] = PY[[k]][idx]
  }
  L_Y_mat[L_Y_mat < 1e-15] = 1e-15
  return((log(L_Y_mat)))
}

get_loglikehood_point_Y3 = function(Y,MissY,b1,b2,b3,B,Z){
  n = dim(Y)[1]
  q = dim(Y)[2]
  XB = Z %*% B 
  cp1 = plogis(XB + matrix(b1, n, q, byrow = TRUE))
  cp2 = plogis(XB + matrix(b2, n, q, byrow = TRUE))
  cp3 = plogis(XB + matrix(b3, n, q, byrow = TRUE))
  PY = list()
  PY[[1]] = cp1
  PY[[2]] = cp2 - cp1
  PY[[3]] = cp3 - cp2
  PY[[4]] = 1 - cp3
  L_Y_mat = matrix(0, n, q)
  for (k in 1:4) {
    idx = (Y == k)
    L_Y_mat[idx] = PY[[k]][idx]
  }
  L_Y_mat[L_Y_mat < 1e-15] = 1e-15
  L_Y_mat[MissY == 1] = 0
  return((log(L_Y_mat)))
}

get_loglikehood_point_Y4 = function(Y,gamma_Y,B,Z){
  n = dim(Y)[1]
  q = dim(Y)[2]  
  theta_Y = Z %*% B + (matrix(1,n,1) %*% matrix(gamma_Y,1,q))
  L_Y = Y * theta_Y - log(1 + exp(theta_Y))
  return((L_Y))
}


log_mean_exp <- function(x) {
  max_x <- max(x)
  return(max_x + log(mean(exp(x - max_x))))
}

get_WAIC_data3 = function(A, Y1, Y2, Y3, Y4, MissY1, MissY3, 
                          b11_mcmc, b12_mcmc, b13_mcmc, b14_mcmc, 
                          b21_mcmc, b22_mcmc, b31_mcmc, b32_mcmc, b33_mcmc, 
                          gamma4_mcmc, B1_mcmc, B2_mcmc, B3_mcmc, B4_mcmc, 
                          alpha_mcmc, Z_mcmc) {
  
  n = dim(A)[1]
  q1 = ncol(Y1); q2 = ncol(Y2); q3 = ncol(Y3); q4 = ncol(Y4)
  k = dim(Z_mcmc)[3]
  len = dim(alpha_mcmc)[1]
  
  L_A = array(NA, c(len, n, n))
  L_Y1 = array(NA, c(len, n, q1))
  L_Y2 = array(NA, c(len, n, q2))
  L_Y3 = array(NA, c(len, n, q3))
  L_Y4 = array(NA, c(len, n, q4))
  
  for (i in 1:len) {
    Zi = matrix(Z_mcmc[i, , ], n, k)
    
    L_A[i,,] = get_loglikehood_point_A(A, alpha_mcmc[i,], Zi)
    L_Y1[i,,] = get_loglikehood_point_Y1(Y1, MissY1, b11_mcmc[i,], b12_mcmc[i,], b13_mcmc[i,], b14_mcmc[i,], matrix(B1_mcmc[i,,], k, q1), Zi)
    L_Y2[i,,] = get_loglikehood_point_Y2(Y2, b21_mcmc[i,], b22_mcmc[i,], matrix(B2_mcmc[i,,], k, q2), Zi)
    L_Y3[i,,] = get_loglikehood_point_Y3(Y3, MissY3, b31_mcmc[i,], b32_mcmc[i,], b33_mcmc[i,], matrix(B3_mcmc[i,,], k, q3), Zi)
    L_Y4[i,,] = get_loglikehood_point_Y4(Y4, gamma4_mcmc[i,], matrix(B4_mcmc[i,,], k, q4), Zi)
  }
  
  calc_lppd <- function(L_array) {
    apply(L_array, c(2, 3), log_mean_exp)
  }
  
  calc_pwaic <- function(L_array) {
    apply(L_array, c(2, 3), var)
  }
  
  lppd_A = calc_lppd(L_A); pwaic_A = calc_pwaic(L_A)
  lppd_Y1 = calc_lppd(L_Y1); pwaic_Y1 = calc_pwaic(L_Y1)
  lppd_Y2 = calc_lppd(L_Y2); pwaic_Y2 = calc_pwaic(L_Y2)
  lppd_Y3 = calc_lppd(L_Y3); pwaic_Y3 = calc_pwaic(L_Y3)
  lppd_Y4 = calc_lppd(L_Y4); pwaic_Y4 = calc_pwaic(L_Y4)
  
  ut = upper.tri(lppd_A, diag = FALSE)
  wa = sum(lppd_A[ut]) - sum(pwaic_A[ut])
  
  sum_valid <- function(lppd, pwaic) {
    valid = !is.infinite(lppd) & !is.na(lppd) & (lppd != 0)
    sum(lppd[valid]) - sum(pwaic[valid])
  }
  
  wy = sum_valid(lppd_Y1, pwaic_Y1) + 
    sum_valid(lppd_Y2, pwaic_Y2) + 
    sum_valid(lppd_Y3, pwaic_Y3) + 
    sum_valid(lppd_Y4, pwaic_Y4)
  
  return(-2 * (wa + wy))
}

get_DIC_data3 = function(A, Y1, Y2, Y3, Y4, MissY1, MissY3, 
                         b11_mcmc, b12_mcmc, b13_mcmc, b14_mcmc, b21_mcmc,
                         b22_mcmc, b31_mcmc, b32_mcmc, b33_mcmc, gamma4_mcmc,
                         B1_mcmc, B2_mcmc, B3_mcmc, B4_mcmc, alpha_mcmc, Z_mcmc) {
  
  n = dim(Y1)[1]
  q1 = dim(Y1)[2]
  q2 = dim(Y2)[2]
  q3 = dim(Y3)[2]
  q4 = dim(Y4)[2]
  k = dim(Z_mcmc)[3]  
  len = dim(alpha_mcmc)[1]
  

  alpha_mean = apply(alpha_mcmc, 2, mean)
  b11_mean = apply(b11_mcmc, 2, mean)
  b12_mean = apply(b12_mcmc, 2, mean)
  b13_mean = apply(b13_mcmc, 2, mean)
  b14_mean = apply(b14_mcmc, 2, mean)
  b21_mean = apply(b21_mcmc, 2, mean)
  b22_mean = apply(b22_mcmc, 2, mean)
  b31_mean = apply(b31_mcmc, 2, mean)
  b32_mean = apply(b32_mcmc, 2, mean)
  b33_mean = apply(b33_mcmc, 2, mean)
  gamma4_mean = apply(gamma4_mcmc, 2, mean)
  Z_mean = apply(Z_mcmc, c(2,3), mean)
  B1_mean = apply(B1_mcmc, c(2,3), mean)
  B2_mean = apply(B2_mcmc, c(2,3), mean)
  B3_mean = apply(B3_mcmc, c(2,3), mean)
  B4_mean = apply(B4_mcmc, c(2,3), mean)
  
  L_AY_mean = get_loglikehood_data3(A, Y1, Y2, Y3, Y4, MissY1, MissY3,
                                    b11_mean, b12_mean, b13_mean, b14_mean, 
                                    b21_mean, b22_mean, b31_mean, b32_mean, b33_mean,
                                    gamma4_mean, B1_mean, B2_mean, B3_mean, B4_mean,
                                    alpha_mean, Z_mean)
  D_hat = -2 * L_AY_mean
  
  bar_D = numeric(len)
  for (i in 1:len) {
    alphai = alpha_mcmc[i, ]
    b11i = b11_mcmc[i, ]
    b12i = b12_mcmc[i, ]
    b13i = b13_mcmc[i, ]
    b14i = b14_mcmc[i, ]
    b21i = b21_mcmc[i, ]
    b22i = b22_mcmc[i, ]
    b31i = b31_mcmc[i, ]
    b32i = b32_mcmc[i, ]
    b33i = b33_mcmc[i, ]
    gamma4i = gamma4_mcmc[i, ]
    Zi = matrix(Z_mcmc[i, , ], nrow = n, ncol = k)
    B1i = matrix(B1_mcmc[i, , ], nrow = k, ncol = q1)
    B2i = matrix(B2_mcmc[i, , ], nrow = k, ncol = q2)
    B3i = matrix(B3_mcmc[i, , ], nrow = k, ncol = q3)
    B4i = matrix(B4_mcmc[i, , ], nrow = k, ncol = q4)
    
    L_AYi = get_loglikehood_data3(A, Y1, Y2, Y3, Y4, MissY1, MissY3,
                                  b11i, b12i, b13i, b14i, b21i, b22i, 
                                  b31i, b32i, b33i, gamma4i,
                                  B1i, B2i, B3i, B4i, alphai, Zi)
    bar_D[i] = -2 * L_AYi
  }
  
  # DIC = 2*mean(D) - D_hat
  return(2 * mean(bar_D) - D_hat)
}

create_folds <- function(n, k = 5) {
  shuffled_indices <- sample(1:n)
  fold_indices <- split(shuffled_indices, ceiling(seq_along(shuffled_indices) / (length(shuffled_indices) / k)))
  folds_list <- list()
  for (i in 1:k) {
    test_indices <- fold_indices[[i]]
    train_indices <- unlist(fold_indices[-i])
    folds_list[[i]] <- list(
      train = train_indices,
      test = test_indices
    )
  }
  return(folds_list)
}

get_loglikehood_data3 = function(A,Y1,Y2,Y3,Y4,MissY1,MissY3,b11,b12,b13,b14,b21,
                                 b22,b31,b32,b33,gamma4,B1,B2,B3,B4,alpha,Z){
  L_A = get_loglikehood_A(A,alpha,Z)
  L_Y1 = get_loglikehood_Y1(Y1,MissY1,b11,b12,b13,b14,B1,Z)
  L_Y2 = get_loglikehood_Y2(Y2,b21,b22,B2,Z)
  L_Y3 = get_loglikehood_Y3(Y3,MissY3,b31,b32,b33,B3,Z)
  L_Y4 = get_loglikehood_Y4(Y4,gamma4,B4,Z)
  return(L_A + L_Y1 + L_Y2 + L_Y3 + L_Y4)
}