initf <- function(){
  list(v=t(matrix(c(0.1,0.9,0.9,1,0.8,0.8,0.1,0.9,0.1,0.1,0.3,0.1,0.7,0.7,0.7,0.6),16,4)),
       mu_vi_wt_tmp=c(0.1,0.1,0.1,0.1,0.3,0.1),
       mu_vi_g=t(matrix(c(0.1,0.1,0.1,0.1,0.3,0.1,0.6),7,3)),
       a=rep(0.1,16),
       b=rep(0.1,11),
       e_effs=rep(0.1,22),

       sigma_n=0.1,
       sigma_n2=0.1,
       sigma_p=0.1,
       r_p_tmp=rep(0,14),
       r_p=t(matrix(0,16,3)),
       met_est=matrix(1.0,2,3)       )
}