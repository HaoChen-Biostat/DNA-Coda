########################### Summary of Support Recovery ############################
#----------------------------------------------------------------------------------#
# ROC_PR: 
#     Computes True Positive Rate (TPR) and False Positive Rate (FPR) of a delta estimator.
# summary_ROC_PR: 
#     Computes TPR and FPR along the solution path (different tuning parameters).
# compute_pr_norm:
#     Computes norm-based errors and average ROC metrics (TPR/FPR) for replicate estimates.
#----------------------------------------------------------------------------------#

ROC_PR <- function(pred, true){
  #-----------------------------------#
  # Input: 
  #       pred, (p * p) matrix, 
  #           the estimator of delta;
  #       true, (p * p) matrix,
  #           true delta. 
  #
  # Output:
  #       TPR, numeric, 
  #            true positive rate;
  #       FPR, numeric, 
  #            false positive rate.
  #-----------------------------------#
  Upper_true <- true[upper.tri(true)]
  Upper_pred <- pred[upper.tri(pred)]
  TP <- sum(Upper_pred[Upper_true != 0] != 0)
  FP <- sum(Upper_pred[Upper_true == 0] != 0)
  FN <- sum(Upper_pred[Upper_true != 0] == 0)
  TN <- sum(Upper_pred[Upper_true == 0] == 0)
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  return(c(TPR, FPR))
}


summary_ROC_PR = function(path, true){
  #---------------------------------------------------------------------#
  # Input: 
  #       path, (nlambda * 1) list, 
  #           the estimators of delta under different tuning parameters;
  #       true, (p * p) matrix, 
  #           true delta. 
  #
  # Output:
  #       result, (nlambda * 2) matrix, 
  #           TPR and FPR under different tuning parameters.
  #---------------------------------------------------------------------# 
  
  result <- data.frame(t(sapply(path, function(pred)ROC_PR(pred, true))))
  colnames(result) <- c('TPR', 'FPR')
  return(result)
}






########################### Summary of Norm and ROC Metrics ###########################
#---------------------------------------------------------------------------------------#
# compute_pr_norm:
#     computes summary statistics (norms and ROC metrics) for a set of estimated matrices,
#     automatically applies thresholding to all estimates along the path and CV-selected estimates.
#     If no threshold is specified, defaults to 0 (i.e., no change).
#
# Input:
#     est_path        : list of length T (number of replicates), each element is a list of length L (number of λ values)
#                       containing (p×p) estimated delta matrices along the tuning path.
#     est_cv          : list of length T of (p×p) cross-validated estimates (optimal λ)
#     true_delta      : (p×p) matrix of true delta.
#     threshold       : numeric (default 0); coefficients with absolute value ≤ threshold are set to zero.
#
# Output:
#     A list with three components:
#       $mean_se    : a 7×2 matrix of mean and standard error across replicates for metrics
#                     rows are: L2 norm, L1 norm, Frobenius norm, Infinity norm, Max norm, TPR, FPR.
#       $roc_tpr    : numeric vector of length L+1, average true positive rates along est_path, last entry = 1.0.
#       $roc_fpr    : numeric vector of length L+1, average false positive rates along est_path, last entry = 1.0.
#---------------------------------------------------------------------------------------#

compute_pr_norm <- function(est_path, est_cv, true_delta, threshold = 0) {
  T <- length(est_path)
  L <- length(est_path[[1]])
  
  # (1) Thresholding all estimates
  for (t in seq_len(T)) {
    for (l in seq_len(L)) {
      m <- est_path[[t]][[l]]
      m[abs(m) <= threshold] <- 0
      est_path[[t]][[l]] <- m
    }
    m_cv <- est_cv[[t]]
    m_cv[abs(m_cv) <= threshold] <- 0
    est_cv[[t]] <- m_cv
  }
  
  # (2) Norm summaries using thresholded CV estimates
  M <- matrix(nrow = 7, ncol = T)
  for (t in seq_len(T)) {
    mat <- est_cv[[t]]
    M[, t] <- c(
      norm(mat - true_delta, "2"),
      norm(mat - true_delta, "1"),
      norm(mat - true_delta, "F"),
      norm(mat - true_delta, "I"),
      norm(mat - true_delta, "M"),
      ROC_PR(mat, true_delta)
    )
  }
  
  #mean_se <- round(cbind(rowMeans(M), apply(M, 1, sd)), 5)
  se_fn <- function(x) sd(x) / sqrt(T)
  mean_se <- round(cbind(
    mean = rowMeans(M),
    se   = apply(M, 1, se_fn)
  ), 5)
  
  rownames(mean_se) <- c("S_norm","L1_norm","F_norm","∞_norm","m_norm","TPR","FPR")
  colnames(mean_se) <- c("mean", "se")
  
  # (3) ROC metrics along full path using thresholded path estimates
  roc_list <- lapply(est_path, function(path) summary_ROC_PR(path, true_delta))
  tpr_mat <- sapply(roc_list, `[[`, "TPR")
  fpr_mat <- sapply(roc_list, `[[`, "FPR")
  
  roc_tpr <- c(rowMeans(tpr_mat), 1.0)
  roc_fpr <- c(rowMeans(fpr_mat), 1.0)
  
  return(list(
    mean_se = mean_se,
    roc_tpr = roc_tpr,
    roc_fpr = roc_fpr
  ))
}

