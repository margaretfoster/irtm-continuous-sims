
mse <- function(ypred, ytrue, aggregate=TRUE, root=FALSE){
  if(length(dim(ypred)) == 3) 
    ypred <- apply(ypred, c(1,2), mean)

  ytrue_sd <- apply(ytrue, 2, sd)
  ypred_sd <- apply(ypred, 2, sd)
  if (aggregate)
    ret <- mean(colMeans((ytrue/ytrue_sd-ypred/ypred_sd)^2))
  else
    ret <- (ytrue/ytrue_sd-ypred/ypred_sd)^2
  if (root)
    ret <- sqrt(ret)

  return(ret)
}

## unstandardized, for CFA (not preds not divided by variance)

mse_unstd <- function(ypred, ytrue, aggregate=TRUE, root=FALSE){
  if(length(dim(ypred)) == 3) 
    ypred <- apply(ypred, c(1,2), mean)
  
  if (aggregate)
    ret <- mean(colMeans((ytrue-ypred)^2))
  else
    ret <- (ytrue-ypred)^2
  if (root)
    ret <- sqrt(ret)
  
  return(ret)
}

sbc <- function(true_theta, theta){
  if(length(dim(theta))==3)
    return(sbc_md(true_theta, theta))
  else
    return(sbc_1d(true_theta, theta))
}

sbc_1d <- function(true_theta, theta){
  tdiff <- nrow(theta) - length(true_theta) + 1
  stat <- rowSums(true_theta < theta[tdiff:nrow(theta), ])
  stat
}

sbc_md <- function(true_theta, theta){
  tdiff <- nrow(theta) - nrow(true_theta) + 1
  true_sum <- rowSums(true_theta)
  theta_sum <- apply(theta[tdiff:nrow(theta),,], c(1,3), sum)    
  stat <- rowSums(true_sum < theta_sum[, ])
  stat
}

sbc_ks <- function(true_theta, theta){
    S <- dim(theta)[length(dim(theta))]
    stat <- sbc(true_theta, theta)
    ks.test(sample(S, nrow(true_theta), replace=T), stat)
}


gew_stat <- function(j){2*pnorm(-abs(geweke.diag(j)$z)) < 0.05}
ess_stat <- function(j){effectiveSize(j)}
rhat_stat <- function(j){rhat(j)}

geweke_fcn <- function(res){
  mean(apply(res, c(1,2), gew_stat), na.rm=T)
}

ess_fcn <- function(res){
  mean(apply(res, c(1,2), ess_stat), na.rm=T)
}

rhat_fcn <- function(res){
  mean(apply(res, c(1,2), rhat_stat), na.rm=T)
}

cred_int <- function(x, a){
  sx <- sort(x)
  hi <- sx[floor((1 - a) * length(x))]
  lo <- sx[floor(a * length(x))]
  c(lo, hi)
}

coverage <- function(res, true){
  cis <- apply(res, c(1,2), cred_int, 0.05)
  mod_res <- mean((true >= cis[1,,]) & (true <= cis[2,,]))
  mod_res
}

standardize_theta <- function(theta, Sigma){
  theta_std <- t(t(theta)/sqrt(diag(Sigma)))
  theta_std
}

theta_lambda_traceplots <- function(irt, i=NULL, k=NULL){
  theta = irt$theta
  lambda = irt$lambda
  if (is.null(i))
    i = sample(1:nrow(theta), 1)
  if (is.null(k))
    k = sample(1:nrow(lambda), 1)
  n = length(theta[i, 1, ])
  D = ncol(theta)
  par(mfcol=c(D, 3))
  for(d in 1:D){
    plot(1:n, theta[i, d, ], xlab = 'sample', ylab=paste('theta', d))
  }
  for(d in 1:D){
    plot(1:n, lambda[k, d, ], xlab = 'sample', ylab=paste('lambda', d))
  }
  for(d in 1:D){
    plot(1:n, theta[i, d, ] * lambda[k, d, ], xlab = 'sample', ylab=paste('theta * lambda', d))
  }
  par(mfrow=c(1,1))
}
