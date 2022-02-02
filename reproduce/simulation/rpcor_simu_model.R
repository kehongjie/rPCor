

# -------------------- generate X --------------------
Sigma_X <- ar1_cor(p, rho1)


## generate MVN/pareto X
if (heavy_tail) {
  X <- matrix(NA, n, p)
  for(j in 1:p){
    X[,j] <- rpareto(n, 1, shape = 1)
  }
  if(rho1!=0) {
    X <- X %*% chol(Sigma_X)
  }
  
} else {
  X <- mvrnormArma(n, rep(0,p), Sigma_X)
}


# -------------------- construct beta matrix --------------------
beta_matrix <- matrix(0, p, q)


if (!heavy_tail) { ## normal case
  
  if (p==5000) { ## p=5000 or p=1000
    num_x <- 100
    num_y <- 50
    num_edge <- 200
  } else { ## else p=1000
    num_x <- 5
    num_y <- 5
    num_edge <- 10
  }
  
  x_pos <- sample(1:p, num_x) ## number of signal X
  y_pos <- sample(1:q, num_y) ## number of signal Y
  
  for (i in x_pos) {
    beta_matrix[i,sample(y_pos,1)] <- 1
  }
  for (j in y_pos) {
    beta_matrix[sample(x_pos,1),j] <- 1
  }
  idx <- sum(beta_matrix)
  while(idx < num_edge) { ## number of signal edges
    one_x <- sample(x_pos, 1)
    one_y <- sample(y_pos, 1)
    if (beta_matrix[one_x, one_y]==0) {
      beta_matrix[one_x, one_y] <- 1
      idx <- idx + 1 
    }
  }
  
}

if (heavy_tail) { ## heavy tailed case
  beta_matrix[c(1, 2, 3), 1] <- 1
  beta_matrix[c(4, 5, 6), 2] <- 1
  beta_matrix[c(7, 8, 9), 3] <- 1
}


beta_matrix <- signal_strength * beta_matrix


# -------------------- generate Y --------------------
if (block==1) {
  ## construct block-wise correlation matrix Sigma_Y
  all_image_block_corr <- list()
  for (h in 1:R) {
    if (h > (R-R_prime)) { ## singular Y that does not belong to any group
      one_block <- matrix(rho2_btw,V,V) ## Y's that do not belong to any group
      diag(one_block) <- 1
    } else {
      if (sig_ytype == "CP") {
        ## for each block: compound symmetric
        one_block <- matrix(rho2,V,V)
        for (j in 1:V) {
          one_block[j,j] <- 1
        }
      }
      if (sig_ytype == "AR1") {
        ## for each block: AR1
        one_block <- matrix(0,V,V)
        for (i in 1:V) {
          for (j in 1:V) {
            one_block[i,j] <- rho2^(abs(i-j))
          }
        }
      }
    }
    ## save correlation matrix for h-th block
    all_image_block_corr[[h]] <- one_block
  }
  Sigma_Y <- as.matrix(bdiag(all_image_block_corr)) 
  
  ## between-group correlation
  if (rho2_btw != 0) {
    full <- matrix(rho2_btw,q,q)
    all_block <- list()
    for (r in 1:R) {
      all_block[[r]] <- matrix(rho2_btw,V,V)
    }
    diagnonal <- as.matrix(bdiag(all_block))
    off_diag <- full - diagnonal
    Sigma_Y <- Sigma_Y + off_diag
  }
} else {
  Sigma_Y <- ar1_cor(q, rho2)
}

Y <- matrix(NA,n,q)
for (i in 1:n) {
  mean_Y <- X[i,] %*% beta_matrix
  one_Y <- mvrnormArma(1, mean_Y, Sigma_Y)
  Y[i,] <- one_Y
}


# -------------------- standardize and index --------------------
X_std <- scale(X)
Y_cen <- scale(Y, scale = F)
Y_std <- scale(Y)

signal_grp_index <- matrix(0, G, R)
for (i in 1:G) {
  for (j in 1:R) {
    sub_beta <- beta_matrix[((i-1)*S+1):((i-1)*S+S), ((j-1)*V+1):((j-1)*V+V)]
    if (sum(sub_beta > 0) > 0) { signal_grp_index[i,j] <- 1}
  }
}

true_group_X <- c(rep(1:(G-G_prime),each=S), rep(0,G_prime*S)) ## the group information for each X
true_group_Y <- c(rep(1:(R-R_prime),each=V), rep(0,R_prime*V)) ## the group information for each Y

true_x_pos <- which(rowSums(beta_matrix)!=0)
true_y_pos <- which(colSums(beta_matrix)!=0)

signal_pos <- which(beta_matrix != 0)
noise_pos <- which(beta_matrix == 0)

