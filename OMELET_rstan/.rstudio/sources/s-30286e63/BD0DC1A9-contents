library(rstan)
library(shinystan)
library(bayesplot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

mkpath <- function(load_dir,fname){
  return(paste(load_dir,fname,sep=""))
}

mkdata <- function(load_dir,fname){
  path_now <- mkpath(load_dir,fname)
  list_tmp <- read.table(path_now,header=FALSE,sep=",")
  return(as.matrix(list_tmp))
}

mkdata_char <- function(load_dir,fname){
  path_now <- mkpath(load_dir,fname)
  list_tmp <- read.table(path_now,header=FALSE,sep="\t")
  return(as.matrix(list_tmp))
}

mkdata_vec <- function(load_dir,fname){
  path_now <- mkpath(load_dir,fname)
  list_tmp <- read.table(path_now,header=FALSE,sep=",")
  return(as.vector(list_tmp$V1))
}

WAIC <- function(fit){
  log_lik <- rstan::extract(fit)$log_lik
  num_r <- dim(log_lik)[2]
  lppd_r <- c()
  p_waic_r <- c()
  waic_r <- c()
  for(r in 1:num_r){
    log_lik_r <- log_lik[,r,]
    p_waic_r[r] <- mean(colMeans(log_lik_r^2) - colMeans(log_lik_r)^2)
    lppd_r[r] <- mean(log(colMeans(exp(log_lik_r))))
    waic_r[r] <- -lppd_r[r]+p_waic_r[r]
  }
  waic <- sum(waic_r)
  return(list(waic=waic,waic_r=waic_r,lppd_r=lppd_r,p_waic_r=p_waic_r))
}

my_stan <- function(save_dir,data_path,initf_path,smodel,
                        iter,warmup,thin,adapt_delta,max_treedepth,
                        c_v,c_e,v_max){
  data <- mkdata_now(data_path,c_v,c_e,v_max)
  save(data,file=paste(save_dir,'/data.RData',sep=""))
  source(initf_path)
  dir.create(save_dir)
  fit = rstan::stan(smodel,
                    data=data,
                    init=initf,
                    chains=4,iter=iter,warmup = warmup,thin=thin,
                    sample_file=paste(save_dir,'/sample',sep=""),
                    diagnostic_file = paste(save_dir,'/diagnostic',sep=""),
                    control = list(adapt_delta = adapt_delta,max_treedepth = max_treedepth))
  file.copy(paste0("./",smodel),
            paste0(save_dir,"/",smodel))
  file.copy(initf_path,
            paste0(save_dir,"/initf.R"))
  waic <- WAIC(fit)
  # launch_shinystan(fit)
  save(fit,file=paste(save_dir,'/fit.RData',sep=""))
  save(waic,file=paste(save_dir,'/waic.RData',sep=""))
  print(waic$waic)
  return(list(fit=fit,waic=waic))
}

my_stan_transform <- function(fit,data,myfun_now){
  ms <- rstan::extract(fit)
  iter <- nrow(ms$a)
  
  source(myfun_now)
  mu_vi_wt_tmp_ <- calc_mu_vi_wt(ms,data)
  
  # mu_vi
  mu_vi <- array(0,dim=c(iter,data$num_g,data$num_ri))
  mu_vi[,1,] <- mu_vi_wt_tmp_
  mu_vi[,2:data$num_g,] <- ms["mu_vi_g"][[1]]
  # mu_v
  mu_v <- array(0,dim=c(iter,data$num_g,data$num_rc))
  for(i in 1:iter){
    tmp <- t(data$Nu %*% t(mu_vi[i,,]))
    mu_v[i,,] <- tmp[,data$idx_calc]
  }
  
  # mean_v
  sum_v <- matrix(0,nrow = iter,ncol=data$num_rc)
  for(g in 1:data$num_g){
    sum_v <- sum_v + mu_v[,g,] * (data$idx_g[g,2]-data$idx_g[g,1]+1)
  }
  mean_v <- sum_v/data$num_smpl
  
  # b
  b_tmp <- ms["b"][[1]]
  idx_b <- rep(1,data$num_rc)
  idx_b[data$is_irrev==1] <- 0
  idx_b[data$is_irrev==1] = 0
  idx_tmp <- which(idx_b==1)
  b <- array(0,dim=c(iter,data$num_rc))
  b[,idx_tmp] <- -b_tmp
  
  # r_p_all
  r_p_all <- array(0,dim=c(iter,data$num_g,data$num_rc))
  if(!is.null(ms$r_p_tmp)){
    r_p <- ms["r_p"][[1]]
    r_p_tmp <- ms["r_p_tmp"][[1]]
    sum_r <- matrix(0,nrow = iter,ncol=2)
    for(g in 1:(data$num_g-1)){
      r_p_all[,g,] = r_p[,g,];
      for(p in 1:data$num_p){
        sum_r[,p] <- sum_r[,p] + r_p[,g,data$idx_p[p]]*sum( data$rna[ data$idx_p[p], data$idx_g[g,1]:data$idx_g[g,2] ] )
      }
    }
    r_p_all[,data$num_g,data$idx_p_] = r_p_tmp;
    for(p in 1:data$num_p){
      r_p_all[,data$num_g,data$idx_p[p]] = -sum_r[,p] / sum( data$rna[data$idx_p[p], data$idx_g[data$num_g,1]:data$idx_g[data$num_g,2]] );
    }
  }else if(is.null(ms$r_p)){
    r_p_all <- ms["r_p"][[1]]
  }
  
  # enz
  enz_pred <- ms$enz_pred
  rna_pred <- array(0,dim=c(iter,data$num_rx,data$num_smpl))
  if(!is.null(ms$rna_pred)){
    rna_pred <- ms$rna_pred
  }

  # met_est
  c_out <- array(0,dim=c(iter,data$num_met_est,data$num_g))
  if(!is.null(ms$met_est)){
    sum_c <- matrix(0,nrow = iter,ncol=data$num_met_est)
    met_est <- ms$met_est
    for(c in 1:data$num_met_est){
      for(g in 1:(data$num_g-1)){
        sum_c[,c] <- sum_c[,c] + met_est[,c,g]*(data$idx_g[g,2]-data$idx_g[g,1]+1);
      }
      for(i in 1:iter){
        c_out[i,c,data$num_g] <- (data$num_smpl-sum_c[i,c]) / (data$idx_g[data$num_g,2]-data$idx_g[data$num_g,1]+1)
      }
    }
  }

  # output
  ms_out <- ms
  ms_out$mu_vi <- mu_vi
  ms_out$mu_v <- mu_v
  ms_out$mean_v <- mean_v
  ms_out$b <- b
  ms_out$r_p_all <- r_p_all
  ms_out$enz_out <- enz_pred
  ms_out$c_out <- c_out
  
  return(ms_out)
}

output_data <- function(ms,save_dir){
  param_list <- names(ms)
  for(par in param_list){
    ms_now <- ms[par][[1]]
    iter <- dim(ms_now)[1]
    ndim <- length(dim(ms_now))
    if(ndim==1){
      write(ms_now,file=paste0(save_dir,"/",par,".txt"),ncolumns = 1)
    }else if(ndim==2){
      write(t(ms_now),file=paste0(save_dir,"/",par,".txt"),ncolumns = dim(ms_now)[2])
    }else{
      for(i in 1:dim(ms_now)[2]){
        ms_now2 <- ms_now[,i,]
        write(t(ms_now2),file=paste0(save_dir,"/",par,".txt"),ncolumns = dim(ms_now2)[2],append = TRUE)
      }
    }
  }
}

dens_flux <- function(pars_plot,ms,data_path,col_list1,col_list2,save_dir){
  for(par in pars_plot){
    ms_now <- ms[par][[1]]
    min_now <- floor(min(ms_now))
    max_now <- ceiling(max(ms_now))
    idx_include <- mkdata_vec(data_path,"/idx_include.txt")
    
    if(par=="v"||par=="a"||par=="b"||par=="r_p_all"){
      rxn_names <- mkdata_vec(data_path,"/rxn_names_include.txt")
      num_sub <- length(rxn_names)
    }else if(par=="mu_vi"){
      rxn_names <- mkdata_vec(data_path,"/rxn_names_indflux.txt")
      num_sub <- length(rxn_names)
    }
    
    pdf(paste0(save_dir, "/dens_",par,".pdf"),width=5,height=num_sub*1.5)
    par(mar = c(2,4,2,2))#lower,left,upper,right
    par(mfrow=c(num_sub,1))
    # if parameters are 3 dimensions
    if(length(dim(ms_now))>2 && num_sub==dim(ms_now)[3]){
      for(r in 1:num_sub){
        ms_now_now <- ms_now[,,r]
        for(g in 1:dim(ms_now)[2]){
          d <- density(ms_now_now[,g])
          if(g==1){
            plot(d,xlim=c(min_now,max_now),main=rxn_names[r],xlab="",cex=5)
            polygon(d,col=col_list2[g],border=col_list1[g]) 
          }else{
            par(new=T)
            polygon(d,col=col_list2[g],border=col_list1[g]) 
            # plot(d,xlim=c(min_now,max_now),cex=5,add=T,axis=F)
          }
        }
      }
    }else if(num_sub==dim(ms_now)[2]){
      for(i in 1:ncol(ms_now)){
        d <- density(ms_now[,i])
        m_now <- median(ms_now[,i])
        sd_now <- sd(ms_now[,i])
        main_now <- paste0(rxn_names[i]," (median=", round(m_now,2), ", sd=", round(sd_now,2),")" )
        plot(d,xlim=c(min_now,max_now),main=main_now,xlab="",cex=5)
        polygon(d,col="grey")
      }
    }
    dev.off()
  }
}


trace_param <- function(fit,save_dir,par,w_now,h_now){
  ms <- rstan::extract(fit)
  # param_list <- names(ms)
  # for(par in param_list){
  # if(par != "log_lik" && par != "y" && par !="x" && par != "y_eff" &&
  # par != "enz_pred" && par != "enz_pred2"
  # ){
  dim_tmp <- dim(ms[[par]])
  ### traceplot
  pdf(paste0(save_dir, "/traceplot_", par ,".pdf"),width=w_now,height=h_now)
  plot_now <- traceplot(fit, pars=par,inc_warmup=T)
  print(plot_now)
  dev.off()
  # }
  # }
}

ac_param <- function(fit,save_dir,par,w_now,h_now){
  ms <- rstan::extract(fit)
  # param_list <- names(ms)
  # for(par in param_list){
  # if(par != "log_lik" && par != "y" && par !="x" && par != "y_eff" &&
  # par != "enz_pred" && par != "enz_pred2"
  # ){
  dim_tmp <- dim(ms[[par]])
  ### traceplot
  pdf(paste0(save_dir, "/acplot_", par ,".pdf"),width=w_now,height=h_now)
  plot_now <- stan_ac(fit, pars=par,separate_chains = F)
  print(plot_now)
  dev.off()
  # }
  # }
}

output_summary <- function(fit,save_dir){
  tmp <- summary(fit)$summary
  write.csv(tmp,paste0(save_dir,'/summary_out.csv'))
}
