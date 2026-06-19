#For simulate 1

analyze_H <- function(res, cn, burn, thin, nrun, k) {
  
  replation <- length(res)
  
  H_net <- array(0, c(replation, length(cn), nrun))
  H_Y   <- array(0, c(replation, length(cn), nrun))
  H_jo  <- array(0, c(replation, length(cn), nrun))
  
  for (i in seq_len(replation)) {
    for (j in seq_along(cn)) {
      H_net[i, j, ] <- res[[i]]$out$net[[j]]$H
      H_Y[i, j, ]   <- res[[i]]$out$Y[[j]]$H
      H_jo[i, j, ]  <- res[[i]]$out$net_Y[[j]]$H
    }
  }
  
  thining <- seq(burn + thin, nrun, thin)
  
  H_net_thin <- H_net[, , thining, drop = FALSE]
  H_Y_thin   <- H_Y[, , thining, drop = FALSE]
  H_jo_thin  <- H_jo[, , thining, drop = FALSE]
  
  H_net_mean <- apply(H_net_thin, c(1, 2), getmode)
  H_Y_mean   <- apply(H_Y_thin, c(1, 2), getmode)
  H_jo_mean  <- apply(H_jo_thin, c(1, 2), getmode)
  
  net_table <- lapply(seq_along(cn), function(j) {
    table(H_net_mean[, j])
  })
  
  Y_table <- lapply(seq_along(cn), function(j) {
    table(H_Y_mean[, j])
  })
  
  jo_table <- lapply(seq_along(cn), function(j) {
    table(H_jo_mean[, j])
  })
  
  names(net_table) <- paste0("cn=", cn)
  names(Y_table)   <- paste0("cn=", cn)
  names(jo_table)  <- paste0("cn=", cn)
  
  return(list(
    
    H_net_mean = H_net_mean,
    H_Y_mean   = H_Y_mean,
    H_jo_mean  = H_jo_mean,
    
    net_rate = get_rate(H_net_mean, k),
    net_bias = get_bias(H_net_mean, k),
    
    Y_rate = get_rate(H_Y_mean, k),
    Y_bias = get_bias(H_Y_mean, k),
    
    jo_rate = get_rate(H_jo_mean, k),
    jo_bias = get_bias(H_jo_mean, k),
    
    net_table = net_table,
    Y_table   = Y_table,
    jo_table  = jo_table
    
  ))
}

analyze_alpha <- function(res, cn) {
  
  replation <- length(res)
  
  ALPHA_BIAS_JO <- vector("list", length(cn))
  ALPHA_BIAS_NET <- vector("list", length(cn))
  
  summary_result <- data.frame(
    cn = cn,
    NET_RMSE_Mean = NA_real_,
    NET_RMSE_SD   = NA_real_,
    JO_RMSE_Mean  = NA_real_,
    JO_RMSE_SD    = NA_real_
  )
  
  for (i in seq_along(cn)) {
    
    n <- cn[i]
    
    alpha_bias_jo  <- matrix(0, replation, n)
    alpha_bias_net <- matrix(0, replation, n)
    
    for (j in seq_len(replation)) {
      
      alpha_true <- res[[j]]$ture[[i]]$alpha
      
      alpha_bias_jo[j, ] <-
        apply(res[[j]]$out$net_Y[[i]]$alpha, 2, mean) -
        alpha_true
      
      alpha_bias_net[j, ] <-
        apply(res[[j]]$out$net[[i]]$alpha, 2, mean) -
        alpha_true
    }
    
    ALPHA_BIAS_JO[[i]] <- alpha_bias_jo
    ALPHA_BIAS_NET[[i]] <- alpha_bias_net
    
    rmse_net <- sqrt(rowMeans(alpha_bias_net^2))
    rmse_jo  <- sqrt(rowMeans(alpha_bias_jo^2))
    
    summary_result$NET_RMSE_Mean[i] <- mean(rmse_net)
    summary_result$NET_RMSE_SD[i]   <- sd(rmse_net)
    
    summary_result$JO_RMSE_Mean[i] <- mean(rmse_jo)
    summary_result$JO_RMSE_SD[i]   <- sd(rmse_jo)
  }
  
  names(ALPHA_BIAS_JO)  <- paste0("cn=", cn)
  names(ALPHA_BIAS_NET) <- paste0("cn=", cn)
  
  return(list(
    ALPHA_BIAS_JO  = ALPHA_BIAS_JO,
    ALPHA_BIAS_NET = ALPHA_BIAS_NET,
    summary = summary_result
  ))
}


analyze_gamma <- function(res, cn, cq) {
  
  replation <- length(res)
  
  GAMMA_BIAS_JO <- vector("list", length(cn))
  GAMMA_BIAS_Y  <- vector("list", length(cn))
  
  summary_result <- data.frame(
    cn = cn,
    cq = cq,
    Y_RMSE_Mean  = NA_real_,
    Y_RMSE_SD    = NA_real_,
    JO_RMSE_Mean = NA_real_,
    JO_RMSE_SD   = NA_real_
  )
  
  for (i in seq_along(cn)) {
    
    q <- cq[i]
    
    gamma_bias_jo <- matrix(0, replation, q)
    gamma_bias_y  <- matrix(0, replation, q)
    
    for (j in seq_len(replation)) {
      
      gamma_true <- res[[j]]$ture[[i]]$gamma
      
      gamma_bias_jo[j, ] <-
        apply(res[[j]]$out$net_Y[[i]]$gamma, 2, mean) -
        gamma_true
      
      gamma_bias_y[j, ] <-
        apply(res[[j]]$out$Y[[i]]$gamma, 2, mean) -
        gamma_true
    }
    
    GAMMA_BIAS_JO[[i]] <- gamma_bias_jo
    GAMMA_BIAS_Y[[i]]  <- gamma_bias_y
    
    rmse_y  <- sqrt(rowMeans(gamma_bias_y^2))
    rmse_jo <- sqrt(rowMeans(gamma_bias_jo^2))
    
    summary_result$Y_RMSE_Mean[i]  <- mean(rmse_y)
    summary_result$Y_RMSE_SD[i]    <- sd(rmse_y)
    
    summary_result$JO_RMSE_Mean[i] <- mean(rmse_jo)
    summary_result$JO_RMSE_SD[i]   <- sd(rmse_jo)
  }
  
  names(GAMMA_BIAS_JO) <- paste0("cn=", cn)
  names(GAMMA_BIAS_Y)  <- paste0("cn=", cn)
  
  return(list(
    GAMMA_BIAS_JO = GAMMA_BIAS_JO,
    GAMMA_BIAS_Y  = GAMMA_BIAS_Y,
    summary       = summary_result
  ))
}


analyze_B <- function(res, cn, cq, thining) {
  
  replation <- length(res)
  
  de_B_Y  <- matrix(0, length(cn), replation)
  de_B_JO <- matrix(0, length(cn), replation)
  
  for (i in seq_along(cn)) {
    
    q <- cq[i]
    
    for (j in seq_len(replation)) {
      
      B_true <- res[[j]]$ture[[i]]$B
      BB_true <- t(B_true) %*% B_true
      
      ## ----- Y -----
      BB_hat_Y <- matrix(0, q, q)
      
      for (mm in seq_along(thining)) {
        
        B_hat_mm <- res[[j]]$out$Y[[i]]$B[[mm]]
        
        if (is.vector(B_hat_mm)) {
          B_hat_mm <- matrix(B_hat_mm, nrow = 1)
          B_hat_mm <- rbind(B_hat_mm, 0)
        }
        
        if (dim(B_hat_mm)[1] == 1) {
          B_hat_mm <- 0 * B_hat_mm
        }
        
        BB_hat_Y <- BB_hat_Y + t(B_hat_mm) %*% B_hat_mm
      }
      
      de_B_Y[i, j] <-
        mean((BB_hat_Y / length(thining) - BB_true)^2)
      
      ## ----- JO -----
      BB_hat_JO <- matrix(0, q, q)
      
      for (mm in seq_along(thining)) {
        
        B_hat_mm <- res[[j]]$out$net_Y[[i]]$B[[mm]]
        
        if (is.vector(B_hat_mm)) {
          B_hat_mm <- matrix(B_hat_mm, nrow = 1)
          B_hat_mm <- rbind(B_hat_mm, 0)
        }
        
        if (dim(B_hat_mm)[1] == 1) {
          B_hat_mm <- 0 * B_hat_mm
        }
        
        BB_hat_JO <- BB_hat_JO + t(B_hat_mm) %*% B_hat_mm
      }
      
      de_B_JO[i, j] <-
        mean((BB_hat_JO / length(thining) - BB_true)^2)
    }
  }
  
  summary_result <- data.frame(
    cn = cn,
    cq = cq,
    
    Y_RMSE_Mean = apply(sqrt(de_B_Y), 1, mean),
    Y_RMSE_SD   = apply(sqrt(de_B_Y), 1, sd),
    
    JO_RMSE_Mean = apply(sqrt(de_B_JO), 1, mean),
    JO_RMSE_SD   = apply(sqrt(de_B_JO), 1, sd)
  )
  
  return(list(
    de_B_Y  = de_B_Y,
    de_B_JO = de_B_JO,
    summary = summary_result
  ))
}


analyze_Z <- function(res, cn, thining) {
  
  replation <- length(res)
  
  de_Z_Y   <- matrix(0, length(cn), replation)
  de_Z_net <- matrix(0, length(cn), replation)
  de_Z_jo  <- matrix(0, length(cn), replation)
  
  for (i in seq_along(cn)) {
    
    n <- cn[i]
    
    for (j in seq_len(replation)) {
      
      Z_true <- res[[j]]$ture[[i]]$Z
      ZZ_true <- Z_true %*% t(Z_true)
      
      ## ----- Y -----
      ZZ_hat_Y <- matrix(0, n, n)
      
      for (mm in seq_along(thining)) {
        Z_hat_mm <- res[[j]]$out$Y[[i]]$Z[[mm]]
        ZZ_hat_Y <- ZZ_hat_Y + Z_hat_mm %*% t(Z_hat_mm)
      }
      
      de_Z_Y[i, j] <-
        mean((ZZ_hat_Y / length(thining) - ZZ_true)^2)
      
      ## ----- NET -----
      ZZ_hat_net <- matrix(0, n, n)
      
      for (mm in seq_along(thining)) {
        Z_hat_mm <- res[[j]]$out$net[[i]]$Z[[mm]]
        ZZ_hat_net <- ZZ_hat_net + Z_hat_mm %*% t(Z_hat_mm)
      }
      
      de_Z_net[i, j] <-
        mean((ZZ_hat_net / length(thining) - ZZ_true)^2)
      
      ## ----- JO -----
      ZZ_hat_jo <- matrix(0, n, n)
      
      for (mm in seq_along(thining)) {
        Z_hat_mm <- res[[j]]$out$net_Y[[i]]$Z[[mm]]
        ZZ_hat_jo <- ZZ_hat_jo + Z_hat_mm %*% t(Z_hat_mm)
      }
      
      de_Z_jo[i, j] <-
        mean((ZZ_hat_jo / length(thining) - ZZ_true)^2)
    }
  }
  
  summary_result <- data.frame(
    cn = cn,
    
    Y_RMSE_Mean   = apply(sqrt(de_Z_Y), 1, mean),
    Y_RMSE_SD     = apply(sqrt(de_Z_Y), 1, sd),
    
    NET_RMSE_Mean = apply(sqrt(de_Z_net), 1, mean),
    NET_RMSE_SD   = apply(sqrt(de_Z_net), 1, sd),
    
    JO_RMSE_Mean  = apply(sqrt(de_Z_jo), 1, mean),
    JO_RMSE_SD    = apply(sqrt(de_Z_jo), 1, sd)
  )
  
  return(list(
    de_Z_Y   = de_Z_Y,
    de_Z_net = de_Z_net,
    de_Z_jo  = de_Z_jo,
    summary  = summary_result
  ))
}