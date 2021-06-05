library(rstan)
# library(shinystan)
# library(bayesplot)
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

mkdata_multi <- function(load_dir,c_v,c_e,v_max,grp,pair){
  int_list_tmp <- mkdata_char(load_dir,"/int_list.txt")
  is_irrev_tmp <- mkdata_vec(load_dir,"/is_irrev.txt")
  idx_include_tmp <- mkdata_vec(load_dir,"/idx_include.txt")
  data = list(num_r=as.integer(int_list_tmp[int_list_tmp[,1]=="num_r",2]),
              num_ri=as.integer(int_list_tmp[int_list_tmp[,1]=="num_ri",2]),
              num_rd=as.integer(int_list_tmp[int_list_tmp[,1]=="num_rd",2]),
              num_rc=as.integer(int_list_tmp[int_list_tmp[,1]=="num_rc",2]),
              num_m=as.integer(int_list_tmp[int_list_tmp[,1]=="num_m",2]),
              num_mc=as.integer(int_list_tmp[int_list_tmp[,1]=="num_rd",2]),
              num_b=as.integer(int_list_tmp[int_list_tmp[,1]=="num_b",2]),
              num_ri_wt=as.integer(int_list_tmp[int_list_tmp[,1]=="num_ri_wt",2]),
              num_smpl=as.integer(int_list_tmp[int_list_tmp[,1]=="num_smpl",2]),
              num_g=as.integer(int_list_tmp[int_list_tmp[,1]=="num_g",2]),
              num_p=as.integer(int_list_tmp[int_list_tmp[,1]=="num_p",2]),
              S=mkdata(load_dir,"/S.txt"),
              N=mkdata(load_dir,"/N.txt"),
              Nu=mkdata(load_dir,"/Nu.txt"),
              Ne=mkdata(load_dir,"/Ne.txt"),
              Sp=mkdata(load_dir,"/Sp.txt"),
              Sm=mkdata(load_dir,"/Sm.txt"),
              idx_calc=mkdata_vec(load_dir,"/idx_calc.txt"),
              idx_g=mkdata(load_dir,"/idx_g.txt"),
              enz=mkdata(load_dir,"/enz.txt"),
              sub1=mkdata(load_dir,"/sub1.txt"),
              pro1=mkdata(load_dir,"/pro1.txt"),
              sub2=mkdata(load_dir,"/sub2.txt"),
              pro2=mkdata(load_dir,"/pro2.txt"),
              is_irrev=is_irrev_tmp[idx_include_tmp],
              c_v=c_v,
              c_e=c_e,
              v_max=v_max,
              met_eff=mkdata(load_dir,"/met_eff.txt"))
  # extract only the target group samples
  idx_g_now <- data$idx_g[grp,]
  data$num_g <- length(grp)
  idx_g <- c()
  for(g in 1:data$num_g){
    idx_g <- append(idx_g,idx_g_now[g,1]:idx_g_now[g,2])
  }
  data$num_smpl <- length(idx_g)
  # normalize data
  enz <- matrix(0,data$num_rc,data$num_smpl)
  sub1 <- matrix(1,data$num_rc,data$num_smpl)
  pro1 <- matrix(1,data$num_rc,data$num_smpl)
  sub2 <- matrix(1,data$num_rc,data$num_smpl)
  pro2 <- matrix(1,data$num_rc,data$num_smpl)
  for(r in 1:data$num_rc){
    if(mean(data$enz[r,idx_g]) != 0){
      enz[r,] <- data$enz[r,idx_g]/mean(data$enz[r,idx_g])
    }
    if(mean(exp(data$sub1[r,idx_g]))-1<1e-1){
      sub1[r,] <- log(exp(data$sub1[r,idx_g])/mean(exp(data$sub1[r,idx_g])))
    }
    if(mean(data$pro1[r,idx_g])-1<1e-1){
      pro1[r,] <- log(exp(data$pro1[r,idx_g])/mean(exp(data$pro1[r,idx_g])))
    }
    if(mean(exp(data$sub2[r,idx_g]))-1<1e-1){
      sub2[r,] <- log(exp(data$sub2[r,idx_g])/mean(exp(data$sub2[r,idx_g])))
    }
    if(mean(data$pro2[r,idx_g])-1<1e-1){
      pro2[r,] <- log(exp(data$pro2[r,idx_g])/mean(exp(data$pro2[r,idx_g])))
    }
    }
  data$enz <- enz
  data$sub1 <- sub1
  data$pro1 <- pro1
  data$sub2 <- sub2
  data$pro2 <- pro2
  
  met_eff <- matrix(0,nrow(data$met_eff),data$num_smpl)
  for(i in 1:nrow(data$met_eff)){
    met_eff[i,] <- log(exp(data$met_eff[i,idx_g])/mean(exp(data$met_eff[i,idx_g])))
    # met_eff[i,] <- exp(data$met_eff[i,idx_g])/mean(exp(data$met_eff[i,idx_g]))
  }
  data$met_eff <- met_eff
  
  idx_tmp <- 1
  idx_g2 <- matrix(0,data$num_g,2)
  for(g in 1:data$num_g){
    idx_g2[g,1] <- idx_tmp;
    idx_g2[g,2] <- idx_tmp + idx_g_now[g,2]-idx_g_now[g,1]
    idx_tmp <- idx_tmp + idx_g_now[g,2]-idx_g_now[g,1] + 1
  }
  data$idx_g <- idx_g2
  
  idx_b <- rep(1,data$num_rc)
  idx_b[data$is_irrev==1] = 0
  data$idx_b <- which(idx_b==1)
  
  return(data)
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

my_stan_multi <- function(save_dir,data_path,initf_path,smodel,
                        iter,warmup,thin,adapt_delta,max_treedepth,
                        c_v,c_e,v_max,grp,pair){
  data <- mkdata_multi(data_path,c_v,c_e,v_max,grp,pair)
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
  # save(fit,file=paste(save_dir,'/fit.RData',sep=""))
  save(waic,file=paste(save_dir,'/waic.RData',sep=""))
  print(waic$waic)
  return(list(fit=fit,waic=waic))
}

