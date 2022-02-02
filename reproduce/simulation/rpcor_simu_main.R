### ### generate seed numbers ### ### 
set.seed(1234)
all_seed <- sample(1:1000000, 1000)
### ### ### ### ### ### ### ### ###


dir_source <- "" ## the directory where all the codes are located
dir_output <- "" ## the directory where the output goes

library(BCSub) ## mvrnormArma()
library(Matrix)
library(foreach)
library(doParallel) ## foreach()
library(parallel)
library(corpcor) ## pcor.shrink()
library(lars)
library(SIS)
library(pcalg) ## pcSelect()
library(igraph)
library(mvtnorm) ## rmvt()
library(robustbase) ## Qn()
library(FarmTest) ## huber.cov()
library(ppcor) ## spcor()
library(energy) ## dcor.test()

# source(paste(dir_source, "/screening_FUN.R", sep="")) ## other screening methods
# source(paste(dir_source, "/pcSelect.R", sep="")) ## for implementing PC-simple, from pcalg R pkg

numCores <- 15 ## number of cores to do parallel computing

# args <- commandArgs(TRUE)
args <- c(1,4,0,1,1,2,0)
## Notes: this "args" is for varying the scenarioes of simulation,
##        which includes hyper-parameters (simu_id, rho_idx,  
##        rho1_btw, sig_type, block, signal_strength, heavy_tail).
##        An example is `args <- c(1,4,0,1,1,2,0)`.
##        More details about each parameter can be found below.


simu_id <- as.numeric(args[1]) ## simulation id, also determine the seed number
print(paste("simu_id =", simu_id))
set.seed(all_seed[simu_id]) ## set seed for this repliaction 

# beta_type <- as.numeric(args[2]) ## beta matrix type 
# print(paste("beta_type=",beta_type,sep=""))

rho_idx <- as.numeric(args[2]) ## within-block correlation indicator
all_rho <- matrix(c(0.2,0.2,0.8,0.2,0.2,0.8,0.8,0.8),4,2,byrow=T)
rho1 <- all_rho[rho_idx,1]
rho2 <- all_rho[rho_idx,2]
print(paste("rho is ", c(rho1,rho2), sep=""))

rho1_btw <- rho2_btw <- as.numeric(args[3]) ## between-block correlation
print(paste("rho1_btw=",rho1_btw,sep="")) 

if (as.numeric(args[4])==1) { ## covariance structure(within block), 1 for CP, 2 for AR1
  sig_xtype <- sig_ytype <- "CP" 
} else {sig_xtype <- sig_ytype <- "AR1"}
print(paste("covariance structure is ",sig_xtype,sep=""))

block <- as.numeric(args[5]) ## 1 for block, 0 for no block (for Y)
print(paste("block =", block==1))

signal_strength <- as.numeric(args[6]) ## signal strength
print(paste("signal =", signal_strength))

if (length(args)>=7) {
  heavy_tail <- ifelse(as.numeric(args[7])==1, TRUE, FALSE) ## 1 for heavy tail, 0 for normal
} else {
  heavy_tail <- FALSE 
}
print(paste("heavy_tail =", heavy_tail))


#################################### FUNCTION #################################### 
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  return(rho^exponent)
}

threshold.cov <- function(sigma, threshold = 0.5) {
  sigma_threshold <- matrix(1,nrow(sigma),ncol(sigma))
  sigma_threshold[(sigma != diag(diag(sigma))) & (abs(sigma) < threshold)] <- 0
  return(sigma_threshold)
}

threshold.pvalue <- function(pvalue, threshold=1e-5) {
  pvalue_threshold <- matrix(0,nrow(pvalue),ncol(pvalue))
  pvalue_threshold[(abs(pvalue) < threshold)] <- 1
  return(pvalue_threshold)
}


## FUNCTION: calculate residuals by "lm" function
lm.partial.cor <- function(x, y, z) {
  res1 <- lm(x~z)$residuals
  res2 <- lm(y~z)$residuals
  return(cor(res1,res2))
}

## FUNCTION: calculate residuals in linear regression by matrix multiplication
cal.res <- function(y,X) {
  A <- solve(t(X) %*% X)
  beta <- A %*% t(X) %*% y
  res <- (y - X %*% beta)
  return(res)
}

## FUNCTION: robust sample correlation 
# t <- c*log(p)
# m <- n/2
# the larger the c, the more truncation
robust.scor <- function(x, y, c=0.5, m=100) {
  f_truc <- function(tau){
    mean(  ( (Z^2 + tau^2 ) - abs( tau^2 - Z^2 ))/2  /(tau^2)  ) -  constant
  }
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  t <- c*log(p)
  
  y_mat <- outer(y,y,'-');
  y_diff2 <- as.numeric(y_mat[lower.tri(y_mat)])
  
  constant <- t/m
  
  Z <- y_diff2*y_diff2/2
  
  tau = uniroot(f_truc, interval=c(1e-10,10^10),extendInt="yes")$root
  vary <- abs(mean(pmin(abs(Z),tau)*sign(Z)))
  
  rscov <- varx <- rep(NA,p)
  
  for(j in 1:p){
    xj <- x[,j]
    x_mat <- outer(xj,xj,'-');
    xj_diff2 <- as.numeric(x_mat[lower.tri(x_mat)])
    
    Z <- xj_diff2*y_diff2/2
    tau = uniroot(f_truc, interval=c(1e-10,10^10),extendInt="yes")$root
    rscov[j] <- mean(pmin(abs(Z),tau)*sign(Z))
    
    Z <- xj_diff2*xj_diff2/2
    tau = uniroot(f_truc, interval=c(1e-10,10^10),extendInt="yes")$root
    varx[j] <- abs(mean(pmin(abs(Z),tau)*sign(Z)))
  }      
  
  rsc <- rscov/(sqrt(varx*vary))
  return(rsc)
}

## FUNCTION: robust sample correlation between X (matrix) and y (vector) by FarmTest pkg
rcor.scr <- function(x, y) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  rscov <- varx <- rep(NA,p)
  
  for(j in 1:p){
    xj <- x[,j]
    cov_mat <- huber.cov(cbind(xj,y))
    if(j==1) {vary <- cov_mat[2,2]}
    varx[j] <- cov_mat[1,1]
    rscov[j] <- cov_mat[1,2]
  }   
  rsc <- rscov/(sqrt(varx*vary))
  return(rsc)
}

## FUNCTION: calculate partial correlation by matrix multiplication
partial.cor <- function(x,y,z) {
  res1 <- cal.res(x,z)
  res2 <- cal.res(y,z)
  pcor <- as.numeric(cor(res1, res2))
  # return(as.numeric(cor(res1,res2)))
  return(list(pcor=pcor, res1=res1, res2=res2))
}



rpar.cor <- function(x, y, z, c=2, m=100) {

  XX <- cbind(x, y, z)
  p <- ncol(XX)
  rcov_matrix <- huber.cov(XX)
  P_matrix <- solve(rcov_matrix)
  rpcor <- -P_matrix[1,2] / sqrt(abs(P_matrix[1,1] * P_matrix[2,2]))

  return(sign(rpcor) * min(abs(rpcor), 1))
}


## FUN: calculate robust partial correlation
rpar.cor.own <- function(x, y, z, c=2, m=100) {
  f_truc <- function(tau){
    mean(  ( (Z^2 + tau^2 ) - abs( tau^2 - Z^2 ))/2  /(tau^2)  ) -  constant
  }
  
  XX <- cbind(x, y, z)
  p <- ncol(XX)
  t <- c*log(p-1) 
  constant <- t/m
  
  rcov_matrix <- matrix(NA, p, p)
  for (i in 1:p) {
    x_mat <- outer(XX[,i], XX[,i], '-');
    x_diff2 <- as.numeric(x_mat[lower.tri(x_mat)])
    Z <- x_diff2*x_diff2/2
    tau = uniroot(f_truc, interval=c(1e-10,10^10),extendInt="yes")$root
    rsvar <- abs(mean(pmin(abs(Z),tau)*sign(Z)))
    rcov_matrix[i,i] <- rsvar
    
    for (j in (i+1):p) {
      y_mat <- outer(XX[,j], XX[,j], '-');
      yj_diff2 <- as.numeric(y_mat[lower.tri(y_mat)])
      Z <- yj_diff2*x_diff2/2
      tau = uniroot(f_truc, interval=c(1e-10,10^10),extendInt="yes")$root
      rscov <- mean(pmin(abs(Z),tau)*sign(Z))
      rcov_matrix[i, j] <- rcov_matrix[j, i] <- rscov
    }
    
  }
  
  P_matrix <- solve(rcov_matrix)
  rpcor <- -P_matrix[1,2] / sqrt(abs(P_matrix[1,1] * P_matrix[2,2]))

  return(sign(rpcor) * min(abs(rpcor), 1))
}


## FUNCTION: calculate robust partial corrleation from robust covariance matrix 
rcov2rpor <- function(X) {
  P_matrix <- solve(X)
  rpcor <- -P_matrix[1,2] / sqrt(P_matrix[1,1] * P_matrix[2,2])
  return(rpcor)
}

## FUNCTION: calculate bi-partial correlation 
bipar.cor <- function(x,y,z1,z2) {
  res1 <- cal.res(x,z1)
  res2 <- cal.res(y,z2)
  return(as.numeric(cor(res1,res2)))
}


## FUNCTION: Fisher's Z transformation and p-value for correlation
fisher.ztrans.cor <- function(r,n) { ## r is correlation
  z  <-  0.5*(log(1+r)-log(1-r))
  pvalue <- pnorm(-abs(z), mean=0, sd=(1/sqrt(n-3)))*2
  return(pvalue)
}


## FUNCTION: Fisher's Z transformation and p-value for partial correlation
fisher.ztrans <- function(r,n,k) { ## r is partial correlation
  z  <-  0.5*(log(1+r)-log(1-r))
  test.stat <- sqrt(n-k-3)*z
  pvalue <- pnorm(-abs(test.stat))*2
  return(pvalue)
}

## FUNCTION: implement CIS with `ppcor` R pkg
cis.ppcor <- function(X, Y, sis_pvalue, group_index_X) {
  n <- nrow(X)
  cis_pvalue <- sis_pvalue
  for (g in 1:max(group_index_X)) { ## group_index_Y fixed at 06/30
    pos_x <- which(group_index_X == g)
    if (length(pos_x)>1 & length(pos_x)<(n-2)) {
      for (j in 1:q) {
        comb_data <- cbind(Y[,j], X[,pos_x])
        # onegrp_pcor <- pcor.shrink(comb_data, lambda=0,verbose=F)[1,-1]
        # cis_pvalue[pos_x,j] <- fisher.ztrans(onegrp_pcor, n, length(pos_x)-1)
        cis_pvalue[pos_x,j] <- spcor(comb_data)$p.value[1,-1]
      }
    }
  }
  return(cis_pvalue)
}


## FUNCTION: implement CIS with the codes from the author
cis.screening <- function(X, Y, sis_pvalue, group_index_X) {
  n <- nrow(X)
  q <- ncol(Y)
  # cov_y <- var(Y) ## added for cis.screening, cov or var?
  cov_y <- rep(NA, q)
  for (j in 1:q) {
    cov_y[j] <- var(Y[,j])
  }
  
  cis_pvalue <- sis_pvalue
  for (g in 1:max(group_index_X)) { ## group_index_Y fixed at 06/30
    pos_x <- which(group_index_X == g)
    if (length(pos_x)>1 & length(pos_x)<(n-2)) {
      for (j in 1:q) {
        comb_data <- cbind(Y[,j], X[,pos_x])
        # onegrp_pcor <- pcor.shrink(comb_data, lambda=0,verbose=F)[1,-1]
        # cis_pvalue[pos_x,j] <- fisher.ztrans(onegrp_pcor, n, length(pos_x)-1)
        aa2 <- pcor.shrink(comb_data, verbose=F, lambda=0)
        aa3 <- pvar.shrink(comb_data, verbose=F, lambda=0, lambda.var=0)
        aa <- aa2[1,-1]
        onegrp_pcor <- aa/(cov_y[j]*sqrt((1-aa*aa2[-1,1])/aa3[1]))
        
        cis_pvalue[pos_x,j] <- fisher.ztrans(onegrp_pcor, n, length(pos_x)-1)
      }
    }
  }
  return(cis_pvalue)
}

## FUNCTION: rPCor algorithm for screening on X&Y with the partial correlation order keeping increasing.
##           Conditioning on both X and Y simultaneously 
##           Selecting conditioning sets based on:
##              correlation threshold (X) and block structure (Y) if block=TRUE
##              correlation threshold if block=FALSE
## (last updated 10/19/2021)

PC.both <- function(X, Y, block=F, group_index_Y=NULL, signal_ind, true_ind, marg.ind=T, 
                    alpha.x=1e-5, alpha.y=1e-5, alpha0=1e-2, cor.X=NULL, cor.Y=NULL, 
                    xcor.thres=0.8, ycor.thres=0.8, max.reach=5, d=20, 
                    ht=F, c=2, m=100, ht.mtd=1) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  sig_rpcor <- NULL
  noi_rpcor <- NULL
  
  if (is.null(cor.X) | is.null(cor.Y)) {
    if (ht) {
      cor.X <- cov2cor(huber.cov(X))
      cor.Y <- cov2cor(huber.cov(Y))
      # huber.cov
    } else {
      cor.X <- cor(X)
      cor.Y <- cor(Y)
    }
  }
  true_pos <- which(true_ind != 0)
  output_ind <- signal_ind
  if (ht) { ## if heavy-tailed, calculate and save robust covariance 
    rcov_matrix <- matrix(NA, p+q, p+q)
  }
  
  if (marg.ind==TRUE) {
    ## marginal screening (correlation)
    for (i in 1:p) { 
      # print(i)
      pos_y <- which(signal_ind[i,] != 0)
      grp_y <- group_index_Y[pos_y]
      
      if (length(pos_y) > 0) {## run marginal test if this X has signal
        cor.Y <- cor(X[,i], Y[,pos_y])
        pvalue_y <- fisher.ztrans.cor(cor.Y, n)
        edge_ind <- as.numeric(pvalue_y<alpha0)
        
        temp_posy <- pos_y[which(edge_ind!=0)]
        temp_grpy <- group_index_Y[pos_y]
        
        output_ind[i,pos_y] <- edge_ind
      }## end if on run marginal test if this X has signal
    }## end loop on X_i
    # print(paste("marginal screening done, left with ",sum(output_ind)," pairs", sep=""))
    print(paste("marginal screening done, left with ",sum(output_ind)," pairs and ", sum(output_ind[true_pos]), " are TP", sep=""))
  }
  
  
  for (l in 1:max.reach) { ## PC-algorithm from high correlation threshold to low threshold

    new_output_ind <- output_ind ## update once every iteration (order)
    
    for (j in 1:q) {
      
      if (block==TRUE) {
        if (group_index_Y[j]!=0) {
          cor_pos_y <- setdiff(which(group_index_Y==group_index_Y[j]), j)
        } else {
          cor_pos_y <- 0
        }
      } else {
        cor_pos_y <- setdiff(which(abs(cor.Y[j,])>ycor.thres), j)
        if (length(cor_pos_y) > d) {
          topd_pos <- order(abs(cor.Y[j,cor_pos_y]), decreasing = T)[1:d]
          cor_pos_y <- cor_pos_y[topd_pos]
        }
      }
      
      for (i in 1:p) {
        if (output_ind[i,j] == 1) {
          act_pos_x <- which(output_ind[,j]==1)
          act_pos_y <- which(output_ind[i,]==1)
          cor_pos_x <- setdiff(which(abs(cor.X[i,])>xcor.thres), i)
          if (length(cor_pos_x) > d) {
            topd_pos <- order(abs(cor.X[i,cor_pos_x]), decreasing = T)[1:d]
            cor_pos_x <- cor_pos_x[topd_pos]
          }
          
          nb_pos_x <- intersect(act_pos_x, cor_pos_x)
          nb_pos_y <- intersect(act_pos_y, cor_pos_y)
          
          if ((length(nb_pos_x)+length(nb_pos_y)) > 0) {
            Z <- cbind(X[,nb_pos_x], Y[,nb_pos_y])
            num_nb <- ncol(Z)
            if (num_nb >= l) {
              all_subset_ind <- combn(1:num_nb, l)
              cl <- 1
              one_edge <- 1
              while(cl<=ncol(all_subset_ind) & one_edge==1) {
                sub_Z <- Z[,all_subset_ind[,cl]] ## condition set
                if (ht) {
                  if (ht.mtd==1) {
                    one_pcor <- rpar.cor.own(X[,i], Y[,j], sub_Z, c=c, m=m)
                  }
                  if (ht.mtd==2) {
                    one_pcor <- rpar.cor(X[,i], Y[,j], sub_Z)
                  }
                } else {
                  one_pcor <- partial.cor(X[,i], Y[,j], sub_Z)$pcor
                }
                pvalue <- fisher.ztrans(one_pcor, n, l)
                if (pvalue > alpha.x) {one_edge <- 0}
                cl <- cl+1
              }
              new_output_ind[i,j] <- one_edge
            } ## end if number of neighbord greater than l
          } ## end if has neighbord
        } ## end if edge active
      } ## end for on 1:p
    } ## end for on 1:q
    
    output_ind <- new_output_ind
    print(paste(l,"-th order screening done, left with ",sum(output_ind)," pairs and ",sum(output_ind[true_pos]), " are TP", sep=""))
  } ## end loop on l
  return(output_ind)
}  



#################################### parameter setting #################################### 

if (heavy_tail) {
  n <- 200 
  p <- 1000 
  q <- 500 
  S <- 50 
  V <- 50 
  (G <- 20) == p/S 
  (R <- 10) == q/V 
  G_prime <- 0.1*G 
  R_prime <- 0.1*R 
} else {
  ## large example
  n <- 500 ## sample szie
  p <- 5000 ## number of epigenetic regulators 
  q <- 2000 ## number of genes 
  S <- 50 ## group size in epigenetic regulators (not used anymore)
  V <- 50 ## group size in genes
  (G <- 100) == p/S ## number of groups in epigenetic regulators (not used anymore)
  (R <- 40) == q/V ## number of groups in genes
  G_prime <- 0.1*G ## number of singular groups in epigenetic regulators (not used anymore)
  R_prime <- 0.1*R ## number of singular groups in genes
}

max_reach <- 5
xcor_thres <- 0.3
ycor_thres <- 0.3
remove_nb <- TRUE
d <- round(n/(log(n)*2))
c0 <- 2
m0 <- 100

all_alpha <- c(1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5)


#################################### replication #################################### 
all_on <- Sys.time() 

source(paste(dir_source, "/model_setting.R", sep="")) ## generate data
print("data generated")

## calculate X-X and Y-Y (robust) correlation
start_rcor <- Sys.time()
if (heavy_tail) {
  cor_X <- cov2cor(huber.cov(X_std))
  cor_Y <- cov2cor(huber.cov(Y_std))
} else {
  cor_X <- cor(X_std)
  cor_Y <- cor(Y_std)
}
end_rcor <- Sys.time()
# print(end_rcor - start_rcor)

## block detection by igraph package
adj_X <- threshold.cov(cor_X, xcor_thres)
adj_Y <- threshold.cov(cor_Y, ycor_thres)
adj_X_sp <- as(adj_X, "sparseMatrix")
adj_Y_sp <- as(adj_Y, "sparseMatrix")
graph_X <- graph_from_adjacency_matrix(adj_X_sp, mode="undirected")
graph_Y <- graph_from_adjacency_matrix(adj_Y_sp, mode="undirected")
group_index_X <-  membership(cluster_walktrap(graph_X))
group_index_Y <-  membership(cluster_walktrap(graph_Y))


num_method <- 6
all_method_rep <- matrix(0, num_method, 7)
colnames(all_method_rep) <- c("modelsize","TP","FP","row.kept","column.kept","row.TP","column.TP")
rownames(all_method_rep) <- c("SIS","PC-simple","CIS","rPCor","rPCor-module","DC")
all_method_alpha <- lapply(vector("list", length(all_alpha)), function(x)all_method_rep) 
## Note: `all_method_alpha` saves all the results. It is a list of length 7, and each element of 
##        the list is the output for one specific alpha value. For each alpha value, the output
##        is in the format of a matrix where each row is a method.  

#################################### SIS #################################### 
sis_start <- Sys.time()
sis_pvalue <- matrix(NA, p, q)
for (j in 1:q) {
  cor_y <- as.vector(cor(X_std, Y[,j]))
  pvalue_y <- fisher.ztrans.cor(cor_y, n)
  sis_pvalue[,j] <- pvalue_y
}
for (k in 1:length(all_alpha)) {
  sis_ind <- threshold.pvalue(sis_pvalue, threshold=all_alpha[k])
  
  all_method_alpha[[k]][1,] <- c(sum(sis_ind!=0), sum(sis_ind[signal_pos]), sum(sis_ind[noise_pos]),
                                 sum(rowSums(sis_ind)>0), sum(colSums(sis_ind)>0),
                                 sum(rowSums(sis_ind)[true_x_pos]>0), sum(colSums(sis_ind)[true_y_pos]>0))
  print(all_method_alpha[[k]][1,])
}
sis_end <- Sys.time()
print("SIS finishes in:")
print(sis_end - sis_start)


#################################### CIS ####################################
# calculate semi-partial correlation (CIS) by ppcor package
cis_start <- Sys.time()
cis_pvalue <- cis.ppcor(X_std, Y, sis_pvalue, group_index_X) ## ppcor pkg
# cis_pvalue <- cis.screening(X_std, Y, sis_pvalue, group_index_X) ## corpcor pkg 

for (k in 1:length(all_alpha)) {
  cis_ind <- threshold.pvalue(cis_pvalue, threshold=all_alpha[k])
  all_method_alpha[[k]][3,] <- c(sum(cis_ind!=0), sum(cis_ind[signal_pos]), sum(cis_ind[noise_pos]),
                                 sum(rowSums(cis_ind)>0), sum(colSums(cis_ind)>0),
                                 sum(rowSums(cis_ind)[true_x_pos]>0), sum(colSums(cis_ind)[true_y_pos]>0))
  print(all_method_alpha[[k]][3,])
}

cis_end <- Sys.time()
print("CIS finishes in")
print(cis_end-cis_start)


#################################### rPCor ####################################
pc_start <- Sys.time()

for (k in 1:length(all_alpha)) {
  alpha_x <- alpha_y <- alpha0 <- all_alpha[k]
  if (heavy_tail) {
    marg_ind <- threshold.pvalue(rsis_pvalue, threshold=all_alpha[k])
  } else {
    marg_ind <- threshold.pvalue(sis_pvalue, threshold=all_alpha[k])
  }
  
  pctop_ind <- PC.both(X=X_std, Y=Y_std, signal_ind=marg_ind, true_ind=beta_matrix, marg.ind=F,
                       alpha.x=alpha_x, alpha.y=alpha_y, alpha0=alpha0, 
                       cor.X=cor_X, cor.Y=cor_Y,  
                       xcor.thres=xcor_thres, ycor.thres=ycor_thres,
                       max.reach=max_reach, d=d, ht=heavy_tail, c=c0, m=m0)
  all_method_alpha[[k]][4,] <- c(sum(pctop_ind!=0), sum(pctop_ind[signal_pos]), sum(pctop_ind[noise_pos]),
                                 sum(rowSums(pctop_ind)>0), sum(colSums(pctop_ind)>0),
                                 sum(rowSums(pctop_ind)[true_x_pos]>0), sum(colSums(pctop_ind)[true_y_pos]>0))
  print(all_method_alpha[[k]][4,])
}

pc_end <- Sys.time()
print("rPCor finishes in")
print(pc_end-pc_start)


#################################### rPCor-module ####################################
pc_start <- Sys.time()
for (k in 1:length(all_alpha)) {
  alpha_x <- alpha_y <- alpha0 <- all_alpha[k]
  if (heavy_tail) {
    marg_ind <- threshold.pvalue(rsis_pvalue, threshold=all_alpha[k])
  } else {
    marg_ind <- threshold.pvalue(sis_pvalue, threshold=all_alpha[k])
  }
  
  pccb_ind <- PC.both(X=X_std, Y=Y_std, block=TRUE, group_index_Y=group_index_Y, signal_ind=marg_ind,
                      true_ind=beta_matrix, marg.ind=F,
                      alpha.x=alpha_x, alpha.y=alpha_y, alpha0=alpha0,
                      cor.X=cor_X, cor.Y=cor_Y,  
                      xcor.thres=xcor_thres, ycor.thres=0,
                      max.reach=max_reach, d=d, ht=heavy_tail, c=c0, m=m0)
  all_method_alpha[[k]][5,] <- c(sum(pccb_ind!=0), sum(pccb_ind[signal_pos]), sum(pccb_ind[noise_pos]),
                                 sum(rowSums(pccb_ind)>0), sum(colSums(pccb_ind)>0),
                                 sum(rowSums(pccb_ind)[true_x_pos]>0), sum(colSums(pccb_ind)[true_y_pos]>0))
  print(all_method_alpha[[k]][5,])
}

pc_end <- Sys.time()
print("rPCor-module finishes in")
print(pc_end-pc_start)


#################################### PC-simple ####################################
# parallel version
pc_start <- Sys.time()
for (k in 1:length(all_alpha)) {
  registerDoParallel(numCores)
  print(paste("alpha =", all_alpha[k]))
  
  pc_col_ind <- foreach (j=1:q, .combine='cbind') %dopar% {
    one_col <- rep(0, p)
    pc <- pcSelect(Y_std[, j], X_std, alpha=all_alpha[k])
    one_col[pc$G] <- 1
    one_col
  }
  all_method_alpha[[k]][2,] <- c(sum(pc_col_ind!=0), sum(pc_col_ind[signal_pos]), sum(pc_col_ind[noise_pos]),
                                 sum(rowSums(pc_col_ind)>0), sum(colSums(pc_col_ind)>0), 
                                 sum(rowSums(pc_col_ind)[true_x_pos]>0), sum(colSums(pc_col_ind)[true_y_pos]>0))
  print(all_method_alpha[[k]][2,])
  stopImplicitCluster()
  
}
pc_end <- Sys.time()
print("PC-simple finishes in:")
print(pc_end-pc_start)


#################################### DC ####################################
## Warning: this takes a lot of time to run
dc_start <- Sys.time()

s2_partx <- rep(NA, p)
s2_party <- rep(NA, q)
for (i in 1:p) {
  x_mat <- outer(X_std[,i], X_std[,i], '-')
  x_diff2 <- as.numeric(x_mat[lower.tri(x_mat)])
  s2_partx[i] <- sum(abs(x_diff2))/(n^2)
}
for (j in 1:q) {
  y_mat <- outer(Y_std[,j], Y_std[,j], '-')
  y_diff2 <- as.numeric(y_mat[lower.tri(y_mat)])
  s2_party[j] <- sum(abs(y_diff2))/(n^2)
}

## parallel version 
registerDoParallel(numCores)
dc_pvalue <- foreach (j=1:q, .combine='cbind') %dopar% {
  one_col <- rep(0, p)
  for (i in 1:p) {
    s2 <- s2_partx[i] * s2_party[j]
    one_col[i] <- (1-pnorm(sqrt(dcov(X_std[,i],Y_std[,j])^2*n/s2)))*2
  }
  one_col
}
stopImplicitCluster()

for (k in 1:length(all_alpha)) {
  dc_ind <- threshold.pvalue(dc_pvalue, threshold=all_alpha[k])
  all_method_alpha[[k]][6,] <- c(sum(dc_ind!=0), sum(sis_ind[signal_pos]), sum(dc_ind[noise_pos]),
                                 sum(rowSums(dc_ind)>0), sum(colSums(dc_ind)>0),
                                 sum(rowSums(dc_ind)[true_x_pos]>0), sum(colSums(dc_ind)[true_y_pos]>0))
  print(all_method_alpha[[k]][6,])
}

dc_end <- Sys.time()
print("DC finishes in:")
print(dc_end - dc_start)




all_off <- Sys.time()
print("Everything finishes in")
print(all_off-all_on)

#################################### save the result #################################### 
save(all_method_alpha,
     file = paste(dir_output,
                  "/rPCor_simu",simu_id,"_rhoidx",rho_idx,"_btwrho",10*rho1_btw,
                  "_sigtype",sig_xtype,"_block",block,"_ss",signal_strength, 
                  "_n",n,"_p",p,"_q",q,
                  "_thres",xcor_thres*100,".RData",sep=""))

