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
              idx_p=mkdata_vec(load_dir,"/idx_p.txt"),
              enz=mkdata(load_dir,"/enz.txt"),
              sub=mkdata(load_dir,"/sub1.txt"),
              pro=mkdata(load_dir,"/pro1.txt"),
              rna=mkdata(load_dir,"/rna.txt"),
              is_irrev=is_irrev_tmp[idx_include_tmp],
              c_v=c_v,
              c_e=c_e,
              v_max=v_max,
              enz_eff=mkdata(load_dir,"/enz_eff.txt"),
              rna_eff=mkdata(load_dir,"/rna_eff.txt"),
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
  sub <- matrix(1,data$num_rc,data$num_smpl)
  pro <- matrix(1,data$num_rc,data$num_smpl)
  rna <- matrix(0,data$num_rc,data$num_smpl)
  for(r in 1:data$num_rc){
    if(mean(data$enz[r,idx_g]) != 0){
      enz[r,] <- data$enz[r,idx_g]/mean(data$enz[r,idx_g])
    }
    if(mean(exp(data$sub[r,idx_g]))-1<1e-1){
      sub[r,] <- log(exp(data$sub[r,idx_g])/mean(exp(data$sub[r,idx_g])))
    }
    if(mean(data$pro[r,idx_g])-1<1e-1){
      pro[r,] <- log(exp(data$pro[r,idx_g])/mean(exp(data$pro[r,idx_g])))
    }
    rna[r,] <- data$rna[r,idx_g]/mean(data$rna[r,idx_g])
    }
  enz_eff <- matrix(0,1,data$num_smpl)
  rna_eff <- matrix(0,1,data$num_smpl)
  enz_eff[1,] <- data$enz_eff[1,idx_g]/mean(data$enz_eff[1,idx_g])
  rna_eff[1,] <- data$rna_eff[1,idx_g]/mean(data$rna_eff[1,idx_g])
  data$enz <- enz
  data$sub <- sub
  data$pro <- pro
  data$rna <- rna
  data$enz_eff <- array(enz_eff,dim=c(1,data$num_smpl))
  data$rna_eff <- array(rna_eff,dim=c(1,data$num_smpl))
  
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
  
  idx_p_ <- seq(1,data$num_rc)
  data$idx_p_ <- idx_p_[(idx_p_!=data$idx_p[1]) & (idx_p_!=data$idx_p[2])]
  idx_b <- rep(1,data$num_rc)
  idx_b[data$is_irrev==1] = 0
  idx_b[12] <- 1
  idx_b[16] <- 0
  data$idx_b <- which(idx_b==1)
  
  if(pair){
    data$num_wt <- idx_g2[1,2]-idx_g2[1,1]+1
    data$num_ob <- idx_g2[2,2]-idx_g2[2,1]+1
    data$enz_wt_tmp <- data$enz[,idx_g2[1,1]:idx_g2[1,2]]
    data$enz_ob_tmp <- data$enz[,idx_g2[2,1]:idx_g2[2,2]]
    data$sub_wt <- data$sub[,idx_g2[1,1]:idx_g2[1,2]]
    data$sub_ob <- data$sub[,idx_g2[2,1]:idx_g2[2,2]]
    data$pro_wt <- data$pro[,idx_g2[1,1]:idx_g2[1,2]]
    data$pro_ob <- data$pro[,idx_g2[2,1]:idx_g2[2,2]]
    data$rna_wt <- data$rna[,idx_g2[1,1]:idx_g2[1,2]]
    data$rna_ob <- data$rna[,idx_g2[2,1]:idx_g2[2,2]]
    data$met_eff_wt <- data$met_eff[,idx_g2[1,1]:idx_g2[1,2]]
    data$met_eff_ob <- data$met_eff[,idx_g2[2,1]:idx_g2[2,2]]
    data$enz_eff_wt <- array(data$enz_eff[,idx_g2[1,1]:idx_g2[1,2]],dim=c(1,data$num_wt))
    data$enz_eff_ob <- array(data$enz_eff[,idx_g2[2,1]:idx_g2[2,2]],dim=c(1,data$num_ob))
    data$rna_eff_wt <- array(data$rna_eff[,idx_g2[1,1]:idx_g2[1,2]],dim=c(1,data$num_wt))
    data$rna_eff_ob <- array(data$rna_eff[,idx_g2[2,1]:idx_g2[2,2]],dim=c(1,data$num_ob))
  }
  
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

my_stan_transform <- function(fit,data_path,grp){
  ms <- rstan::extract(fit)
  data <- mkdata_multi(data_path,0.1,0.01,10,grp,FALSE)
  
  iter <- nrow(ms$a)
  
  # mu_vi
  mu_vi <- array(0,dim=c(iter,data$num_g,data$num_ri))
  mu_vi_tmp <- ms["mu_vi_wt_tmp"][[1]]
  v_Pcx <- ms["v_Pcx"][[1]]
  v_Cs <- ms["v_Cs"][[1]]
  mu_vi_g <- ms["mu_vi_g"][[1]]
  mu_vi[,2:data$num_g,] <- mu_vi_g
  if(dim(mu_vi_tmp)[2]==5){
    mu_vi[,1,c(1:4,7)] <- mu_vi_tmp 
  }else{
    mu_vi[,1,1:4] <- mu_vi_tmp
    for(i in 1:iter){
      mu_vi[i,1,7] <- 1-sum(mu_vi_tmp[i,]) 
    }
  }
  
  mu_vi[,1,5] <- v_Pcx
  mu_vi[,1,6] <- v_Cs
  
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
  idx_b[12] <- 1
  idx_b[16] <- 0
  idx_tmp <- which(idx_b==1)
  b <- array(0,dim=c(iter,data$num_rc))
  b[,idx_tmp] <- -b_tmp
  b[,12] <- -b[,12]
  
  # r_p_all
  r_p <- ms["r_p"][[1]]
  r_p_tmp <- ms["r_p_tmp"][[1]]
  r_p_all <- array(0,dim=c(iter,data$num_g,data$num_rc))
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
  
  # enz
  enz <- data$enz
  enz_eff <- data$enz_eff
  enz_out <- array(0,dim=c(iter,data$num_rc,data$num_smpl))
  y <- ms["y"][[1]]
  for(i in 1:iter){
    enz_out[i,13,] <- enz[13,] * enz_eff
    enz_out[i,data$idx_p_,] <- enz[data$idx_p_,]
    ###############
    for(g in 1:data$num_g){
      for(p in 1:data$num_p){
        enz_out[i,data$idx_p[p],data$idx_g[g,1]:data$idx_g[g,2]] <-
          data$rna[data$idx_p[p],data$idx_g[g,1]:data$idx_g[g,2]]*(1+r_p_all[i,g,data$idx_p[p]])
      }
    }
    
  }
  
  # c_OAA, c_Pyr
  sum_c <- matrix(0,nrow = iter,ncol=2)
  c_OAA <- ms$c_OAA
  c_Pyr <- ms$c_Pyr
  c_OAA_out <- array(0,dim=c(iter,data$num_g))
  c_Pyr_out <- array(0,dim=c(iter,data$num_g))
  for(g in 1:(data$num_g-1)){
    sum_c[,1] <- sum_c[,1] + c_OAA[,g]*(data$idx_g[g,2]-data$idx_g[g,1]+1);
    sum_c[,2] <- sum_c[,2] + c_Pyr[,g]*(data$idx_g[g,2]-data$idx_g[g,1]+1);
    c_OAA_out[,g] <- c_OAA[,g]
    c_Pyr_out[,g] <- c_Pyr[,g]
  }
  for(i in 1:iter){
    c_OAA_out[i,data$num_g] <- (data$num_smpl-sum_c[i,1]) / (data$idx_g[data$num_g,2]-data$idx_g[data$num_g,1]+1)
    c_Pyr_out[i,data$num_g] <- (data$num_smpl-sum_c[i,2]) / (data$idx_g[data$num_g,2]-data$idx_g[data$num_g,1]+1)
  }
  
  # output
  ms_out <- ms
  ms_out$mu_vi <- mu_vi
  ms_out$mu_v <- mu_v
  ms_out$mean_v <- mean_v
  ms_out$b <- b
  ms_out$r_p_all <- r_p_all
  ms_out$enz_out <- enz_out
  ms_out$c_OAA_out <- c_OAA_out
  ms_out$c_Pyr_out <- c_Pyr_out
  
  return(ms_out)
}
