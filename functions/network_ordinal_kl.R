library(pgdraw)
library(LaplacesDemon)
library(MASS)

safe_max <- function(x) {
  if (all(is.na(x))) return(-3) # 如果全是NA，返回预设的下限边界
  return(max(x, na.rm = TRUE))
}

safe_min <- function(x) {
  if (all(is.na(x))) return(3)  # 如果全是NA，返回预设的上限边界
  return(min(x, na.rm = TRUE))
}

network_ordinal_kf = function(A,Y,H,nrun,burn,thin){
  n = dim(Y)[1]
  q = dim(Y)[2]
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
  theta_inv = rep(1,H)
  theta_inv_B = rep(1,H)
  B = matrix(rnorm(H * q),nrow = H, ncol = q)
  if(H > 1){
    for (i in 2:H) {
      for (j in 1:(i-1)) {
        B[i,j] = 0
      }
    } 
  }
  Phi = matrix(1,n,q)#文章中的Z
  W = matrix(1,n,q)
  b1<-rtruncnorm(q, a=-3, b=-1.5, mean = 0, sd = 1)
  b2<-rtruncnorm(q, a=-1.5, b=0, mean = 0, sd = 1)
  b3<-rtruncnorm(q, a=0, b=1.5, mean = 0, sd = 1)
  b4<-rtruncnorm(q, a=1.5, b=3, mean = 0, sd = 1)
  
  N_sample = ceiling((nrun - burn)/thin)
  alpha_hat = matrix(0,N_sample,n)
  b1_hat = matrix(0,N_sample,q)
  b2_hat = matrix(0,N_sample,q)
  b3_hat = matrix(0,N_sample,q)
  b4_hat = matrix(0,N_sample,q)
  B_hat = array(0,dim = c(N_sample,H,q))
  Z_hat = array(0,dim = c(N_sample,n,H))
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
      Z_i = Z[-i, , drop = FALSE]
      D_Ai = (D_A[i,])[-i]
      alp_cons = alpha[i] + alpha[-i]
      kappa_Ai = (A[i,])[-i] - 0.5
      sigma_Zi = chol2inv(chol(diag(theta_inv, nrow = H) + t(Z_i) %*% diag(D_Ai,nrow = (n -1)) %*% Z_i + B%*%diag(W_i,nrow = q)%*%t(B)))
      u_Zi = sigma_Zi%*%((t(Z_i) %*% (kappa_Ai - diag(D_Ai,nrow = (n -1)) %*% alp_cons)) + (B%*%diag(W_i,nrow = q)%*%phi_i))
      Z[i,] = mvrnorm(1,u_Zi,sigma_Zi)
    }
    
    #sample b1234
    pb1<-ifelse(Y==1,Phi,NA)
    pb2<-ifelse(Y==2,Phi,NA)
    pb3<-ifelse(Y==3,Phi,NA)
    pb4<-ifelse(Y==4,Phi,NA)
    pb5<-ifelse(Y==5,Phi,NA)
    pb11 <- apply(pb1, 2, safe_max)
    pb20 <- apply(pb2, 2, safe_min)
    pb21 <- apply(pb2, 2, safe_max)
    pb30 <- apply(pb3, 2, safe_min)
    pb31 <- apply(pb3, 2, safe_max)
    pb40 <- apply(pb4, 2, safe_min)
    pb41 <- apply(pb4, 2, safe_max)
    pb50 <- apply(pb5, 2, safe_min)
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
    
    if((run > burn) &((run-burn) %% thin == 0)){
      b1_hat[m,] = b1
      b2_hat[m,] = b2
      b3_hat[m,] = b3
      b4_hat[m,] = b4
      alpha_hat[m,] = alpha
      B_hat[m,,] = B
      Z_hat[m,,] = Z
      m = m + 1
    }    
  }
  output = list("Z" = Z_hat,"B" = B_hat,"alpha" = alpha_hat,"b1" = b1_hat,
                "b2" = b2_hat,"b3" = b3_hat,"b4" = b4_hat)
  return(output)    
}


network_ordinal_kf_getZ = function(A,Y,b1,b2,b3,b4,B,nrun,burn,thin){
  n = dim(A)[1]
  q = dim(Y)[2]
  H = dim(B)[1]
  W = matrix(1,n,q)
  Phi = matrix(1,n,q)
  theta_inv = rep(1/100,H)
  alpha = rnorm(n)
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)
  N_sample = ceiling((nrun - burn)/thin)
  alpha_hat = matrix(0,N_sample,n)
  Z_hat = array(0,dim = c(N_sample,n,H))
  m = 1 
  for (run in 1:nrun){
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
    Theta_w = Z %*% B - Phi
    W = matrix(pgdraw(2,Theta_w),n,q) 
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
    if((run > burn) &((run-burn) %% thin == 0)){
      alpha_hat[m,] = alpha
      Z_hat[m,,] = Z
      m = m+1
    }
  }
  output = list("alpha" = alpha_hat,"Z" = Z_hat)
  return(output)  
}