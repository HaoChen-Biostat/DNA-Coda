rm(list=ls())

library(MASS)
library(huge)
library(mvtnorm)
library(lava)
library(igraph)
library(parallel)
library(foreach)
library(doParallel)
library(iterators)

set.seed(12345)

#dyn.load("dpm.so")#Linux
dyn.load("dpm.dll")# Windows
source("DNA_Coda.R")
source("Summary_roc_norm.R")



n=100
p =20
times<-100
Nlambda=50

#---------------------------------generate differential network(select one):hub/band/block/power-law ------------------------------------------------------------------
#------------------------------ (hub graph) -------------------------------#
G1 <- huge.generator(
  n = n,
  d = p,
  graph = "hub",
  g = 5,
  prob = NULL
)

Theta1 <- G1$theta
Omega_X_0 = Theta1 * sample(c(-1,1), p * p, replace = TRUE) * runif(p * p, 0.5,0.7)
Omega_X_0[lower.tri(Omega_X_0, diag = FALSE)] <- 0
Omega_X_0<-as.matrix(Omega_X_0)
Omega_X_0 <- Omega_X_0 + t(Omega_X_0)
diag(Omega_X_0) = abs(min(eigen(Omega_X_0)$values)) + 0.5
Sigma_X_0 = solve(Omega_X_0)
Omega_X_0 = solve(Sigma_X_0)
Omega_X_0[abs(Omega_X_0) < 10 ^ -4] <- 0
Omega_Y_0 <- Omega_X_0
Omega_Y_0[1:(p/5*1),1:(p/5*1)] <- -1*Omega_Y_0[1:(p/5*1),1:(p/5*1)]
diag(Omega_Y_0) <- diag(Omega_X_0)
Sigma_Y_0 <- solve(Omega_Y_0)
delta_0<- Omega_Y_0-Omega_X_0
delta_0<- as.matrix(delta_0)
#------------------------------ (band graph) -------------------------------#
band_graph = huge.generator(n, p, graph = "band", g = 2, verbose = FALSE)
band_am = as.matrix(band_graph$theta)
band_graph1 = huge.generator(n, p, graph = "band", g = 1, verbose = FALSE)
band_am1 = as.matrix(band_graph1$theta)
Omega_0 = band_am * 0.5 + band_am1 * 0.2+ diag(runif(p, 0.5, 1))
Omega_0 = (abs(min(eigen(Omega_0)$values)) + 1) * diag(p) + Omega_0
Omega_X_0<-Omega_0
Sigma_X_0 = solve(Omega_X_0)
Omega_Y_0 <- Omega_X_0
Omega_Y_0[1:(p/5*1),1:(p/5*1)] <- -1*Omega_Y_0[1:(p/5*1),1:(p/5*1)]
diag(Omega_Y_0) <- diag(Omega_X_0)
Sigma_Y_0 <- solve(Omega_Y_0)
delta_0<- Omega_Y_0-Omega_X_0
delta_0<- as.matrix(delta_0)
#------------------------------ (block graph) -------------------------------#
nblock = 5 
prob = 20 / p 
block_graph = huge.generator(n, p, graph = "cluster",  g = nblock, prob = prob, verbose = FALSE)
block_am = as.matrix(block_graph$theta)
nedge = sum(block_am) / 2; temp = block_am
edge_binary = rbinom(nedge, 1, 0.5)
edge_binary[edge_binary == 1] = 0.7; edge_binary[edge_binary == 0] = 0.5
temp[upper.tri(temp)] = 0; temp[temp != 0] = edge_binary
Omega_0 = temp + t(temp) + diag(runif(p, 0.5,1))
Omega_0 = (abs(min(eigen(Omega_0)$values)) +1) * diag(p) + Omega_0
Omega_X_0<-Omega_0
Omega_Y_0 <- Omega_X_0
Omega_Y_0[1:(p/5*1),1:(p/5*1)] <- -1*Omega_Y_0[1:(p/5*1),1:(p/5*1)]
diag(Omega_Y_0) <- diag(Omega_X_0)
Sigma_Y_0 <- solve(Omega_Y_0)
Sigma_X_0 = solve(Omega_X_0)
delta_0<- Omega_Y_0-Omega_X_0
delta_0<- as.matrix(delta_0)
#------------------------------ (Power-law graph) -------------------------------#
source("powerlaw data generation.R")#This code needs to be modified for different p: each row is scaled by 2(p = 20),3 (p = 40), 4 (p = 60 or 80), or 5(p = 100) respectively
omega_XY_delta<-list()
omega_XY_delta<-generate_networks(p=20)
Omega_X_0<-omega_XY_delta$a1
Omega_Y_0<-omega_XY_delta$a2
delta_0<-omega_XY_delta$diff_network
Sigma_X_0<-solve(Omega_X_0)
Sigma_Y_0<-solve(Omega_Y_0)


#----------------------------------generate compositional data----------------------------------------------------
delta_0_hat_cv_all<-list()
delta_0_hat_all_nlambda<-list()
fit.cv_all<-list()
delta_0_hat_cv_all_oracle<-list()
delta_0_hat_all_nlambda_oracle<-list()
fit.cv_all_oracle<-list()

T_X_all<-list()
T_Y_all<-list()
Z_X_all<-list()
Z_Y_all<-list()

for (i in 1:times) {
  set.seed(i)
  T_X = rmvnorm(n, rep(0, p), as.matrix(Sigma_X_0))
  T_Y = rmvnorm(n, rep(0, p), as.matrix(Sigma_Y_0))
  T_X_all[[i]]<-T_X
  T_Y_all[[i]]<-T_Y
  W_X = exp(T_X)
  W_Y = exp(T_Y)
  X = W_X / rowSums(W_X)
  Y = W_Y / rowSums(W_Y)
  G = diag(p) - matrix(1, p, p) / p
  Z_X = log(X) %*% G
  Z_Y = log(Y) %*% G
  Z_X_all[[i]]<-Z_X
  Z_Y_all[[i]]<-Z_Y
  
}

#---------------------------DNA_Coda_and_Oracle_parallel-------------------------------------------------------------
cl <- makeCluster(10)
registerDoParallel(cl)

start_time <- Sys.time()

result<-foreach(i = 1:times, .errorhandling = 'pass', .packages = c("MASS", "huge","lava", "mvtnorm")) %dopar% {
  set.seed(i)
  #dyn.load("dpm.so")
  dyn.load("dpm.dll")
  Z_X<-Z_X_all[[i]]
  Z_Y<-Z_Y_all[[i]]
  T_X<-T_X_all[[i]]
  T_Y<-T_Y_all[[i]]
  fit.cv_parallel <- DNA_Coda(Z_X,Z_Y,nlambda=Nlambda,tuning="cv",folds=5)
  fit.cv_parallel_oracle<- DNA_Coda(T_X,T_Y,nlambda=50,tuning="cv",folds=5)
  return(list(fit.cv_parallel,fit.cv_parallel_oracle))

}

parallel::stopCluster(cl)
end_time <- Sys.time()
cat("计算时间：", end_time - start_time, "\n")

#------------------------------(DNA_Coda_and_Oracle_result)----------------------------------

DNA_Coda_diff_result_cv = NULL
DNA_Coda_diff_result_cv_oracle= NULL
for (i in 1:times) {
  fit.cv_all[[i]]<-result[[i]][[1]]
  fit.cv_all_oracle[[i]]<-result[[i]][[2]]
}

for (j in 1:times) {
  #DNA_Coda_path_cv_estimator
  fit.cv<-fit.cv_all[[j]]
  delta_0_hat_all_nlambda[[j]]<-fit.cv[["DNA_Coda"]]#path
  delta_0_hat_cv_all[[j]]<-fit.cv[["DNA_Coda"]][[fit.cv[["opt"]]]]#cv
  
  #oracle_path_cv_estimator
  fit.cv_oracle<-fit.cv_all_oracle[[j]]
  delta_0_hat_all_nlambda_oracle[[j]]<-fit.cv_oracle[["DNA_Coda"]]
  delta_0_hat_cv_all_oracle[[j]]<-fit.cv_oracle[["DNA_Coda"]][[fit.cv_oracle[["opt"]]]]
  
  cat(j,"Going on !\n")
  
}

# DNA_Coda_no_threshold
DNA_Coda_res0 <- compute_pr_norm(
  est_path   = delta_0_hat_all_nlambda,
  est_cv     = delta_0_hat_cv_all,
  true_delta = delta_0,
  threshold  = 0
)
DNA_Coda_diff_result_cv <- DNA_Coda_res0$mean_se
DNA_Coda_diff_TPR       <- DNA_Coda_res0$roc_tpr
DNA_Coda_diff_FPR       <- DNA_Coda_res0$roc_fpr
plot(DNA_Coda_diff_FPR, DNA_Coda_diff_TPR, type="l", lwd=2,
     xlab="1-Spec", ylab="Recall", main="DNA_Coda ROC (no threshold)")

# Oracle_no_threshold
Oracle_res0 <- compute_pr_norm(
  est_path   = delta_0_hat_all_nlambda_oracle,
  est_cv     = delta_0_hat_cv_all_oracle,
  true_delta = delta_0,
  threshold  = 0
)
DNA_Coda_diff_result_cv_oracle <- Oracle_res0$mean_se
DNA_Coda_diff_TPR_oracle       <- Oracle_res0$roc_tpr
DNA_Coda_diff_FPR_oracle       <- Oracle_res0$roc_fpr
plot(DNA_Coda_diff_FPR_oracle, DNA_Coda_diff_TPR_oracle, type="l", lwd=2,
     xlab="1-Spec", ylab="Recall", main="Oracle ROC (no threshold)")

# DNA_Coda_threshold
th <- 1e-4
DNA_Coda_res1 <- compute_pr_norm(
  est_path   = delta_0_hat_all_nlambda,
  est_cv     = delta_0_hat_cv_all,
  true_delta = delta_0,
  threshold  = th
)
DNA_Coda_diff_mean_sd_threshold  <- DNA_Coda_res1$mean_se
DNA_Coda_diff_TPR_threshold      <- DNA_Coda_res1$roc_tpr
DNA_Coda_diff_FPR_threshold      <- DNA_Coda_res1$roc_fpr
plot(DNA_Coda_diff_FPR_threshold, DNA_Coda_diff_TPR_threshold, type="l", lwd=2,
     xlab="1-Spec", ylab="Recall", main=paste0("DNA_Coda ROC (threshold=", th, ")"))

# Oracle_threshold
Oracle_res1 <- compute_pr_norm(
  est_path   = delta_0_hat_all_nlambda_oracle,
  est_cv     = delta_0_hat_cv_all_oracle,
  true_delta = delta_0,
  threshold  = th
)
DNA_Coda_diff_mean_sd_threshold_oracle <- Oracle_res1$mean_se
DNA_Coda_diff_TPR_threshold_oracle       <- Oracle_res1$roc_tpr
DNA_Coda_diff_FPR_threshold_oracle       <- Oracle_res1$roc_fpr
plot(DNA_Coda_diff_FPR_threshold_oracle, DNA_Coda_diff_TPR_threshold_oracle, type="l", lwd=2,
     xlab="1-Spec", ylab="Recall", main=paste0("Oracle ROC (threshold=", th, ")"))


DNA_Coda_diff_result_cv
DNA_Coda_diff_mean_sd_threshold
DNA_Coda_diff_result_cv_oracle
DNA_Coda_diff_mean_sd_threshold_oracle

#------------------save_result----------------------------------
save.image("output_20.RData")

