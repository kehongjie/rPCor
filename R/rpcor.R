##' rPCor for variable screening
##'
##' This function will run the \code{rPCor} to perform variable screening
##' on high-dimensional predcitor and high-dimensional response.
##'
##' @import FarmTest
##' @import BCSub
##'
##' @param X: The n*p predictor matrix.
##' @param Y: The n*q response matrix.
##' @param block: Logical; is there any block/module structure in the response
##' variables? If \code{TRUE}, rPCor-module will be run. If \code{FALSE},
##' rPCor-correlation will be run. Default is \code{FALSE} and rPCor-correlation
##' will be run.
##' @param group.index.Y: Group information for response variables (Y). Only need to
##' be specified if \code{block=TRUE}.
##' @param alpha: The alpha threshold. Default is \code{1e-5}.
##' @param xcor.thres: The threshold for between-predictor (robust) correlation.
##' Default is 0.5.
##' @param ycor.thres: The threshold for between-response (robust) correlation. Only
##' need to be specified if \code{block=FALSE}. Default is 0.5.
##' @param max.reach: The maximum order for partial corrlation. Default is 5.
##' @param d: The maximum number of neighbors to be considered as candidates for
##' conditioning set. That is, for predictors \eqn{X_i}, only top d predictors that
##' have the largest marginal (robust) correlation with \eqn{X_i} will be considered
##' for \eqn{X_j}'s conditioning sets. Same for responses if \code{block=FALSE}.
##' @param ht: Logical; is the data heavy-tailed? If \code{TRUE}, robust partial
##' correlation will be calculated. Otherwise, regular partial correlation will be
##' used. Default is \code{FALSE}. Note that \code{ht=TRUE} will significantly
##' increase the computation burden.
##' @param c: Truncation parameter, which will further be used to find the
##' robustification parameter \code{gamma}. See the paper for details. The larger
##' the c, the more truncation. Default is 0.5.
##' @param verbose: Logical; should the results for each partial corrlation order
##' be printed out? Default is \code{TRUE}.
##'
##' @return A p*q 0-1 indicator matrix. Each row is a predictor and each column is
##' a response. For each element of the matrix, 1 means this predictor-response
##' pair is kept after screening while 0 means it's removed.
##'
##' @export
##'
##' @examples
##' \dontrun{
##' set.seed(123)
##'
##' ## FUNCTION: use fisher's transformation to get the p-value for correlation
##' fisher.ztrans.cor <- function(r,n) { ## r is correlation
##'   z  <-  0.5*(log(1+r)-log(1-r))
##'   pvalue <- pnorm(-abs(z), mean=0, sd=(1/sqrt(n-3)))*2
##'   return(pvalue)
##' }
##' 
##' ## FUNCTION: threshold the p-value matrix to get the indicator
##' threshold.pvalue <- function(pvalue, threshold=1e-5) {
##'   indicator <- matrix(0,nrow(pvalue),ncol(pvalue))
##'   indicator[(abs(pvalue) < threshold)] <- 1
##'   return(indicator)
##' }
##' 
##' n <- 200 ## sample szie (1000)
##' p <- 1000 ## number of predictors
##' q <- 500 ## number of responses
##' alpha <- 1e-3 ## alpha threshold
##' 
##' ## simulate data
##' Sigma_X <- toeplitz(0.8^c(0:(p-1)))
##' Sigma_Y <- toeplitz(0.8^c(0:(q-1)))
##' X <- mvrnormArma(n, rep(0,p), Sigma_X)
##' beta_matrix <- matrix(0, p, q)
##' beta_matrix[c(1, 2, 3), 1] <- 1
##' beta_matrix[c(4, 5, 6), 2] <- 1
##' beta_matrix[c(7, 8, 9), 3] <- 1
##' Y <- matrix(NA, n, q)
##' for (i in 1:n) {
##'   mean_Y <- X[i,] %*% beta_matrix
##'   Y[i,] <- mvrnormArma(1, mean_Y, Sigma_Y)
##' }
##' signal_pos <- which(beta_matrix != 0)
##' noise_pos <- which(beta_matrix == 0)
##' X_std <- scale(X)
##' Y_cen <- scale(Y, scale = F)
##' Y_std <- scale(Y)
##' 
##' ## run SIS (marginal screening)
##' sis_start <- Sys.time()
##' sis_pvalue <- matrix(NA, p, q)
##' for (j in 1:q) {
##'   cor_y <- as.vector(cor(X_std, Y[,j]))
##'   pvalue_y <- fisher.ztrans.cor(cor_y, n)
##'   sis_pvalue[,j] <- pvalue_y
##' }
##' sis_ind <- threshold.pvalue(sis_pvalue, threshold=alpha)
##' sis_end <- Sys.time()
##' print("SIS finishes in:")
##' print(sis_end - sis_start)
##' print(paste("Left with ", sum(sis_ind!=0), " pairs and ", 
##'             sum(sis_ind[signal_pos]), " of them are True Positives", sep=""))
##' 
##' ## run rPCor-corr 
##' pc_start <- Sys.time()
##' pc_ind <- rpcor(X=X_std, Y=Y_std, alpha=alpha, xcor.thres=0.5, ycor.thres=0.5)
##' pc_end <- Sys.time()
##' print("PC-thres finishes in")
##' print(pc_end-pc_start)
##' print(paste("Left with ", sum(pc_ind!=0), " pairs and ", 
##'             sum(pc_ind[signal_pos]), " of them are True Positives", sep=""))
##' 
##' }


rpcor <- function(X, Y, block=FALSE, group.index.Y=NULL,
                  alpha=1e-5,
                  xcor.thres=0.5, ycor.thres=0.5, max.reach=5, d=20,
                  ht=FALSE, c=0.5, verbose=TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  signal_ind <- matrix(1, p, q)

  sig_rpcor <- NULL
  noi_rpcor <- NULL

  if (ht) {
    cor.X <- cov2cor(huber.cov(X))
    cor.Y <- cov2cor(huber.cov(Y))
  } else {
    cor.X <- cor(X)
    cor.Y <- cor(Y)
  }

  # true_pos <- which(true_ind != 0)
  output_ind <- signal_ind
  if (ht) { ## if heavy-tailed, calculate and save robust covariance
    rcov_matrix <- matrix(NA, p+q, p+q)
  }

  ## marginal screening (correlation)
  for (i in 1:p) {
    # print(i)
    pos_y <- which(signal_ind[i,] != 0)
    grp_y <- group.index.Y[pos_y]

    if (length(pos_y) > 0) {## run marginal test if this X has signal
      cor.Y <- cor(X[,i], Y[,pos_y])
      pvalue_y <- fisher.ztrans.cor(cor.Y, n)
      edge_ind <- as.numeric(pvalue_y<alpha)

      temp_posy <- pos_y[which(edge_ind!=0)]
      temp_grpy <- group.index.Y[pos_y]

      output_ind[i,pos_y] <- edge_ind
    }## end if on run marginal test if this X has signal
  }## end loop on X_i

  if (verbose) {
    print(paste("Marginal screening done. Left with ", sum(output_ind)," pairs, ",
                sum(rowSums(output_ind)>0), " rows and ",
                sum(colSums(output_ind)>0), " columns", sep=""))
  }

  for (l in 1:max.reach) {
    new_output_ind <- output_ind ## update once every iteration (order)

    for (j in 1:q) {

      if (block==TRUE) {
        if (group.index.Y[j]!=0) {
          cor_pos_y <- setdiff(which(group.index.Y==group.index.Y[j]), j)
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
                  one_pcor <- rpar.cor.own(X[,i], Y[,j], sub_Z, c=c)
                } else {
                  one_pcor <- partial.cor(X[,i], Y[,j], sub_Z)$pcor
                }
                # pvalue <- pcor.ttest(one_pcor, n, l)

                pvalue <- fisher.ztrans(one_pcor, n, l)
                if (pvalue > alpha) {one_edge <- 0}
                cl <- cl+1
              }
              new_output_ind[i,j] <- one_edge
            } ## end if number of neighbord greater than l
          } ## end if has neighbord
        } ## end if edge active
      } ## end for on 1:p
    } ## end for on 1:q

    output_ind <- new_output_ind
    if (verbose) {
      print(paste(l, "-th order screening done. Left with ", sum(output_ind)," pairs, ",
                  sum(rowSums(output_ind)>0), " rows and ",
                  sum(colSums(output_ind)>0), " columns", sep=""))
    }
  } ## end loop on l

  return(output_ind)
  # return(list(sig_rpcor=sig_rpcor, noi_rpcor=noi_rpcor))
}




## FUNCTION: calculate residuals in linear regression by matrix multiplication
cal.res <- function(y,X) {
  A <- solve(t(X) %*% X)
  beta <- A %*% t(X) %*% y
  res <- (y - X %*% beta)
  return(res)
}

## FUNCTION: calculate partial correlation by matrix multiplication
partial.cor <- function(x,y,z) {
  res1 <- cal.res(x,z)
  res2 <- cal.res(y,z)
  pcor <- as.numeric(cor(res1, res2))
  # return(as.numeric(cor(res1,res2)))
  return(list(pcor=pcor, res1=res1, res2=res2))
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


## FUNCTION: robust sample correlation (not used)
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


## FUNCTION: robust partial correlation
rpar.cor.own <- function(x, y, z, c=0.5) {
  f_truc <- function(tau){
    mean(  ( (Z^2 + tau^2 ) - abs( tau^2 - Z^2 ))/2  /(tau^2)  ) -  constant
  }

  XX <- cbind(x, y, z)
  p <- ncol(XX)
  n <- nrow(XX)
  m <- floor(n/2)
  t <- c*log(p-1) ## what is p here?
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









