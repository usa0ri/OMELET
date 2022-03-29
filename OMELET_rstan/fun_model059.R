calc_mu_vi_wt <- function(ms,data){
  mu_vi_wt_tmp <- ms$mu_vi_wt_tmp
  iter <- dim(ms$a)[1]
  mu_vi_wt <- matrix(0,data$num_ri,iter)
  for(i in 1:iter){
    mu_vi_wt[1:6,i] <- mu_vi_wt_tmp[i,1:6];
    mu_vi_wt[7,i] <- 1-data$Nu[data$idx_fixed,] %*% mu_vi_wt[,i];
  }
  return(mu_vi_wt)
}