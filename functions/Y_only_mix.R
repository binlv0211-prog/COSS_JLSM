library(pgdraw)
library(MASS)
library(truncnorm)
library(LaplacesDemon)

Y_mix = function(Y1,Y2,Y3,nrun,burn,thin,delta_n,alpha_H,a_sig,b_sig,a_theta_B,b_theta_B,
                 theta_inf,start_adapt,Hmax,alpha0,alpha1){
  n = dim(Y1)[1]
  q1 = dim(Y1)[2]
  q2 = dim(Y2)[2]
  q3 = dim(Y3)[2]
  q = q1 + q2 + q3
  
  u<-runif(nrun)
  H = Hmax + 1
  Hstar = Hmax
  Z = matrix(rnorm(n * H),nrow = n,ncol = H)  
  theta_inv_B = rep(1,H)
  
  w = rep(1,H)
  zta = rep(1,H)  
  B1 = matrix(rnorm(H * q1),nrow = H, ncol = q1)
  
  B2 = matrix(rnorm(H * q2),nrow = H, ncol = q2)
  
  B3 = matrix(rnorm(H * q3),nrow = H, ncol = q3)
  
  B = cbind(B1,B2,B3)
  
  gamma1 = rep(1,q1)
  gamma2 = rep(1,q2)
  
  Phi = matrix(1,n,q3)#
  W = matrix(1,n,q3)
  b1<-rtruncnorm(q3, a=-3, b=-1.5, mean = 0, sd = 1)
  b2<-rtruncnorm(q3, a=-1.5, b=0, mean = 0, sd = 1)
  b3<-rtruncnorm(q3, a=0, b=1.5, mean = 0, sd = 1)
  b4<-rtruncnorm(q3, a=1.5, b=3, mean = 0, sd = 1)
  
  inv_sigma = rep(1,q1)
  
  N_sample = ceiling((nrun - burn)/thin)
  H_hat = rep(NA,nrun)
  
  Z_hat = list()
  B1_hat = list()
  B2_hat = list()
  B3_hat = list()
  gamma1_hat = matrix(0,N_sample,q1)
  gamma2_hat = matrix(0,N_sample,q2)
  sigma_hat = matrix(0,N_sample,q1)
  b1_hat = matrix(0,N_sample,q3)
  b2_hat = matrix(0,N_sample,q3)
  b3_hat = matrix(0,N_sample,q3)
  b4_hat = matrix(0,N_sample,q3)
  m = 1  
  
  for (run in 1:nrun){
    theta_Y2 = Z %*% B2 + (matrix(1,n,1) %*% matrix(gamma2,1,q2))
    D_Y2 = matrix(pgdraw(1,theta_Y2),nrow = n,ncol = q2)
    
    #sample Phi
    Theta_Y3 = Z %*% B3
    x0 = matrix(-10,n,q3)
    b1_temp = outer(rep(1,n),b1)
    b2_temp = outer(rep(1,n),b2)
    b3_temp = outer(rep(1,n),b3)
    b4_temp = outer(rep(1,n),b4)
    x1 = matrix(10,n,q3)
    phi1<-matrix(rtruncnorm(1, a=x0, b=b1_temp, mean = Theta_Y3, sd = sqrt(1/W)),n,q3)#1
    phi1 = ifelse(Y3==1,phi1,0)
    phi2<-matrix(rtruncnorm(1, a=b1_temp, b=b2_temp, mean = Theta_Y3, sd = sqrt(1/W)),n,q3)#2
    phi2 = ifelse(Y3==2,phi2,0)
    phi3<-matrix(rtruncnorm(1, a=b2_temp, b=b3_temp, mean = Theta_Y3, sd = sqrt(1/W)),n,q3)#3
    phi3 = ifelse(Y3==3,phi3,0)
    phi4<-matrix(rtruncnorm(1, a=b3_temp, b=b4_temp, mean = Theta_Y3, sd = sqrt(1/W)),n,q3)#4
    phi4 = ifelse(Y3==4,phi4,0)
    phi5<-matrix(rtruncnorm(1, a=b4_temp, b=x1, mean = Theta_Y3, sd = sqrt(1/W)),n,q3)#5
    phi5 = ifelse(Y3==5,phi5,0)
    Phi = phi1 + phi2 + phi3 + phi4 + phi5
    #sample W
    Theta_w = Z %*% B3 - Phi
    W = matrix(pgdraw(2,Theta_w),n,q3)     
    
    for(i in 1:n){
      D_Y2i = D_Y2[i,]
      kappa_Y2i = Y2[i,] - 0.5
      
      W_i = W[i,]
      phi_i = Phi[i,]
      
      Sigma_Zi = chol2inv(chol(diag(1, nrow = H) 
                               + B1 %*% diag(inv_sigma,nrow = q1) %*% t(B1)
                               + B2 %*% diag(D_Y2i,nrow = q2) %*% t(B2)
                               + B3 %*% diag(W_i,nrow = q3) %*% t(B3)))
      u_Zi = Sigma_Zi %*% ((B1 %*% diag(inv_sigma,nrow = q1) %*%(Y1[i,] - gamma1))
                           + (B2 %*%(kappa_Y2i - diag(D_Y2i,nrow = q2) %*% gamma2))
                           + (B3 %*%diag(W_i,nrow = q3)%*%phi_i))
      Z[i,] = mvrnorm(1,u_Zi,Sigma_Zi)
    }
    
    
    
    for(j in 1:q1){
      sigma_B1j = chol2inv(chol(diag(theta_inv_B, nrow = H) + inv_sigma[j] * t(Z) %*% Z))
      u_B1j = inv_sigma[j] * sigma_B1j %*% t(Z) %*% (Y1[,j] - gamma1[j])
      B1[,j] = mvrnorm(1, u_B1j, sigma_B1j)    
    }    
    
    for(j in 1:q2){
      D_Y2j = D_Y2[,j]
      sigma_B2j = chol2inv(chol(diag(theta_inv_B, nrow = H) + t(Z) %*% diag(D_Y2j, nrow = n) %*% Z))
      u_B2j = sigma_B2j %*% t(Z) %*% (Y2[,j] - 0.5 - gamma2[j] * D_Y2j)
      B2[,j] = mvrnorm(1, u_B2j, sigma_B2j)
    }
    
    for (j in 1:q3) {
      W_j = W[,j]
      phi_j = Phi[,j]
      sigma_B3j = chol2inv(chol(t(Z)%*%diag(W_j,nrow = n)%*%Z + diag(theta_inv_B, nrow = H)))
      u_B3j = sigma_B3j%*%t(Z)%*%diag(W_j,nrow = n)%*%phi_j
      B3[,j] = mvrnorm(1,u_B3j,sigma_B3j)
    }    
    
    B = cbind(B1,B2,B3)
    
    lhd_spike<-rep(0,H)
    lhd_slab<-rep(0,H)
    for(h in 1:H){
      lhd_spike[h] = exp(sum(log(dnorm(B[h,], mean = 0, sd = theta_inf^(1/2), log = FALSE))))
      lhd_slab[h] = dmvt(x = B[h,], mu=rep(0,q), S=(b_theta_B/a_theta_B)*diag(q), df=2*a_theta_B)
      prob_h = w*c(rep(lhd_spike[h],h),rep(lhd_slab[h],H - h))
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
        v[h] = rbeta(1, shape1 = (Hmax + 1)**(delta_n) + sum(zta == h), shape2 = 1 + sum(zta > h))
      }else{
        v[h] = rbeta(1, shape1 = alpha_H + sum(zta == h), shape2 = 1 + sum(zta > h))
      }
    }
    v[H] = 1
    w[1] = v[1]
    for (h in 2:H){
      w[h] = v[h]*prod(1-v[1:(h-1)])  
    }
    #update inv_theta_B
    for (h in 1:H){
      if (zta[h] <= h){
        theta_inv_B[h] = theta_inf^(-1)
      }
      else{
        theta_inv_B[h] = rgamma(n=1,shape = a_theta_B + 0.5 * q,rate=b_theta_B + 0.5 * t(B[h,]) %*% B[h,])
      }
    }
    #update inv_sigma
    diff_Y1 = apply((Y1 - (matrix(1,n,1) %*% matrix(gamma1,1,q1)) - (Z %*% B1))**2, 2, sum)
    for (j in 1:q1){
      inv_sigma[j] = rgamma(n=1,shape=a_sig + 0.5 * n,rate=b_sig + 0.5 * diff_Y1[j])
    }
    
    si_g1 = 1 / (n * inv_sigma + 1/100)
    sigma_gamma1 = diag(si_g1)#
    #u_gamma_temp = Y - Z %*% B
    u_gamma1 = inv_sigma * si_g1 * (apply(Y1 - Z %*% B1, 2, sum))
    gamma1 = mvrnorm(1, u_gamma1, sigma_gamma1)    
    
    si_g2 = apply(D_Y2, 2, sum) + 1/100
    sigma_gamma2 = diag(1 / si_g2)
    u_gamma2_temp = Y2 - 0.5 - D_Y2 * (Z %*% B2)
    u_gamma2 = sigma_gamma2 %*% (apply(u_gamma2_temp, 2, sum))
    gamma2 = mvrnorm(1, u_gamma2, sigma_gamma2)    
    
    pb1<-ifelse(Y3==1,Phi,NA)
    pb2<-ifelse(Y3==2,Phi,NA)
    pb3<-ifelse(Y3==3,Phi,NA)
    pb4<-ifelse(Y3==4,Phi,NA)
    pb5<-ifelse(Y3==5,Phi,NA)
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
    b10 = runif(q3,pb11,pb20)
    b20 = runif(q3,pb21,pb30)
    b30 = runif(q3,pb31,pb40)
    b40 = runif(q3,pb41,pb50)
    b1 = rtruncnorm(1, a=pb11, b=pb20, mean = b10, sd = 10^5)
    b2 = rtruncnorm(1, a=pb21, b=pb30, mean = b20, sd = 10^5)
    b3 = rtruncnorm(1, a=pb31, b=pb40, mean = b30, sd = 10^5)
    b4 = rtruncnorm(1, a=pb41, b=pb50, mean = b40, sd = 10^5)
    
    diff_Y1 = apply((Y1 - (matrix(1,n,1) %*% matrix(gamma1,1,q1)) - (Z %*% B1))**2, 2, sum)
    for (j in 1:q1){
      inv_sigma[j] = rgamma(n=1,shape=a_sig + 0.5 * n,rate=b_sig + 0.5 * diff_Y1[j])
    }
    
    #print(inv_sigma)
    #update H[t]
    active = which(zta > c(1:H))
    Hstar = length(active)
    if (run >= start_adapt & u[run] <= exp(alpha0 + alpha1 * run)){
      if (Hstar < H - 1){
        # set truncation to Hstar[t] and subset all variables, keeping only active columns
        H = Hstar + 1
        theta_inv_B = c(theta_inv_B[active],theta_inf^(-1))
        w = c(w[active],1-sum(w[active]))
        Z = cbind(Z[,active],rnorm(n,mean=0,sd=sqrt(theta_inf)))
        B1 = rbind(B1[active,],rnorm(q1,mean=0,sd=sqrt(theta_inf)))
        B2 = rbind(B2[active,],rnorm(q2,mean=0,sd=sqrt(theta_inf)))
        B3 = rbind(B3[active,],rnorm(q3,mean=0,sd=sqrt(theta_inf)))
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
        theta_inv_B = c(theta_inv_B,theta_inf^(-1))
        Z = cbind(Z,rnorm(n,mean=0,sd=sqrt(theta_inf)))
        B1 = rbind(B1,rnorm(q1,mean=0,sd=sqrt(theta_inf)))
        B2 = rbind(B2,rnorm(q1,mean=0,sd=sqrt(theta_inf)))
        B3 = rbind(B3,rnorm(q1,mean=0,sd=sqrt(theta_inf)))
        #theta_inv_B = c(theta_inv_B,rgamma(1,a_theta_B,b_theta_B))
        #B = rbind(B,rnorm(q,0,sqrt(theta_inv_B[H])))
        zta = c(zta,H-1)
      }
    }
    
    H_hat[run] = Hstar
    if((run > burn) &((run-burn) %% thin == 0)){
      gamma1_hat[m,] = gamma1
      gamma2_hat[m,] = gamma2
      b1_hat[m,] = b1
      b2_hat[m,] = b2
      b3_hat[m,] = b3
      b4_hat[m,] = b4
      if(Hstar>0){
        B1_hat[[m]] = B1[1:Hstar,, drop=FALSE]
        B2_hat[[m]] = B2[1:Hstar,, drop=FALSE]
        B3_hat[[m]] = B3[1:Hstar,, drop=FALSE]
        Z_hat[[m]] = Z[,1:Hstar, drop=FALSE]
      }else{
        B1_hat[[m]] = B1
        B2_hat[[m]] = B2
        B3_hat[[m]] = B3
        Z_hat[[m]] = Z
      }
      sigma_hat[m,] = 1 / inv_sigma
      m = m+1
    }    
  }
  output = list("H" = H_hat,"B1" = B1_hat,"B2" = B2_hat,"B3" = B3_hat,
                "b1" = b1_hat,"b2" = b2_hat,"b3" = b3_hat,"b4" = b4_hat,
                "gamma1" = gamma1_hat,"gamma2" = gamma2_hat,"Z" = Z_hat,"sigma" = sigma_hat)
  return(output)
}