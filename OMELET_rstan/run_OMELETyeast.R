
source("my_rstan_OMELETyeast.R")

save_dir = "./result/result_yeast"
dir.create("./result")
data_path = "./input_OMELETyeast"
initf_path <- "./initf_OMELETyeast.R"
smodel = "OMELETyeast.stan"
dir.create(save_dir)

# thin <- 2
# adapt_delta <- 0.99
# max_treedepth <- 20
thin <- 1
adapt_delta <- 0.8
max_treedepth <- 10
c_v <- 0.1
c_e <- 0.01
v_max <- 5
grp <- c(1,2,3,4,5)
pair <- FALSE
# res <- my_stan_multi(save_dir,data_path,initf_path,smodel,
#                      20000,17500,thin,adapt_delta,max_treedepth,
#                      c_v,c_e,v_max,
#                      grp,pair)
res <- my_stan_multi(save_dir,data_path,initf_path,smodel,
                     1000,500,thin,adapt_delta,max_treedepth,
                     c_v,c_e,v_max,
                     grp,pair)

##################
ms <- rstan::extract(res$fit)
save(ms,file=paste(save_dir,'/ms.RData',sep=""))
fit <- res$fit
save(fit,file=paste(save_dir,'/fit.RData',sep=""))

col_list <- c("#646464","#005AFF","#FF3200","#03AF7A","#F6AA00")

source("my_plot.R")
pars_plot <- c("v")
for(par in pars_plot){
  dens_flux(par,ms,data_path,col_list,col_list,save_dir)
}
output_data(ms,save_dir)

trace_param(fit,save_dir)
output_summary(fit,save_dir)

# Fig. S3B
plot_violin_flux(ms,col_list,data_path,save_dir)

flux_pca(data_path,data,ms,save_dir)

# Fig. S3C
plot_correlation(ms,col_list,data_path,save_dir)