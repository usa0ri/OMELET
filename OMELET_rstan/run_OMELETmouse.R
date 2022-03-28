
source("my_rstan_OMELETmouse.R")
data_path = "./input_OMELETmouse"
source(paste0(data_path,"/test_OMELETmodel_mkdata.R"))

save_dir = "./result/result_mouse"
dir.create("./result")
dir.create(save_dir)

initf_path <- "./initf_OMELETmouse.R"
source(initf_path)
smodel = paste0(data_path,"/test_OMELETmodel.stan")

# test
thin <- 1
adapt_delta <- 0.8
max_treedepth <- 10
c_v <- 0.1
c_e <- 0.01
v_max <- 5
data <- mkdata_now(data_path,c_v,c_e,v_max)
res <- my_stan(save_dir,data_path,initf_path,smodel,
               1000,500,thin,adapt_delta,max_treedepth,
               c_v,c_e,v_max)
# Settings in Uematsu et al.
#thin <- 2
#adapt_delta <- 0.99
#max_treedepth <- 20
#res <- my_stan_multi(save_dir,data_path,initf_path,smodel,
#                     20000,10000,thin,adapt_delta,max_treedepth,
#                     c_v,c_e,v_max)

# if you want to view the result by shinystan, run the following.
# source("my_rstan_opt.R")
# launch_shinystan(res$fit)

##################
col_list1 <- c("#0063ff","#4dc4ff","#ff0000","#ff8082")
col_list2 <- c("#00000000","#00000000","#00000000","#00000000")

fit <- res$fit
ms <- my_stan_transform(fit,data_path,grp)
save(ms,file=paste(save_dir,'/ms.RData',sep=""))

source("my_plot.R")
pars_plot <- c("v","a","mu_vi","b","r_p_all")
for(par in pars_plot){
  dens_flux(par,ms,data_path,col_list1,col_list2,save_dir)
}
pars_plot <- c("e_allo","e_cofactor")
fname <- "elasticity"
dens_param(pars_plot,ms,data_path,fname,save_dir)

pars_plot <- c("c_OAA_out","c_Pyr_out")
fname <- "c_OAA_Pyr"
dens_param(pars_plot,ms,data_path,fname,save_dir)

output_data(ms,save_dir)
output_summary(fit,save_dir)
trace_param(fit,save_dir)

out <- matrix(0,4*16,10)
for(i in 1:4){
  for(j in 1:16){
    out[16*(i-1)+j,] <- summary(fit)$summary[ paste0("v[",i,",",j,"]"),]
  }
}
write.csv(t(out),paste0(save_dir, "/v_summary.csv"))

# Fig. 4A
dens_flux_all(ms,data_path,col_list1,col_list2,save_dir)
dens_flux_prior_all(ms,data_path,col_list1,col_list2,save_dir)
