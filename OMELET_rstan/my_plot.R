check_sto <- function(ms,data_path,col_list,save_dir){
  is_fluxi <- as.logical(mkdata_vec(data_path,"/is_indflux.txt"))
  is_include <- as.logical(mkdata_vec(data_path,"/is_include.txt"))
  idx_fluxi <- mkdata_vec(data_path,"/idx_indflux.txt")
  idx_include <- mkdata_vec(data_path,"/idx_include.txt")
  
  v <- ms$v
  vi <- ms$mu_vi
  
  Nu <- mkdata(data_path,"/Nu.txt")
  idx_calc <- mkdata_vec(data_path,"/idx_calc.txt")
  N <- Nu[idx_calc,]
  S <- mkdata(data_path,"/S_.txt")
  N_full <- mkdata(data_path,"/N.txt")
  int_list_tmp <- mkdata_char(data_path,"/int_list.txt")
  num_rd <- as.integer(int_list_tmp[int_list_tmp[,1]=="num_rd",2])
  met_names <- mkdata_vec(data_path,"/met_names_int.txt")
  num_sub <- ceiling(sqrt(num_rd))
  
  for(g in 1:num_g){
    vi_now <- vi[,g,]
    v_now <- v[,g,]
    vi_now[,idx_fluxi %in% idx_include] <- v_now[,is_fluxi[idx_include]]
    dev_now <-N%*%t(vi_now)-t(v_now)
    dev_now_rxn <-rowMeans(dev_now)
    
    v_now_new <- N_full%*%t(vi_now)
    v_now_new[idx_include,] <- t(v_now)
    
    Sv_now <- S%*%v_now_new
    Sv_now_m <- rowMeans(Sv_now)
    
    ##################### output
    pdf(paste0(save_dir, "/hist_Sv_group",g,".pdf"),width=10,height=10)
    par(mfrow=c(num_sub,num_sub)) 
    for(m in 1:num_rd){
      Sv_now_now <- Sv_now[m,]
      hist(Sv_now_now,
           col = paste0(col_list[g],"40"),
           border = col_list[g],
           main=met_names[m],
           xlim = c(-1,1),
           xlab = "")
    }
    dev.off()
    
    ##################### normalize by flux size
    dxdt_now <- abs(S)%*%v_now_new
    Sv_now2 <- Sv_now/dxdt_now
    
    ##################### output
    pdf(paste0(save_dir, "/hist_Sv_wt_group",g,".pdf"),width=10,height=10)
    par(mfrow=c(num_sub,num_sub)) 
    for(m in 1:num_rd){
      Sv_wt_now <- Sv_now2[m,]
      hist(Sv_wt_now,
           col = paste0(col_list[g],"40"),
           border = col_list[g],
           main=met_names[m],
           xlim = c(-1,1),
           xlab = "")
    }
    dev.off()
  }
}

check_enz_fitting <- function(ms,grp,data_path,col_list,save_dir){
  int_list_tmp <- mkdata_char(data_path,"/int_list.txt")
  num_rc <- as.integer(int_list_tmp[int_list_tmp[,1]=="num_rc",2])
  rxn_names <- mkdata_vec(data_path,"/rxn_names_include.txt")
  
  enz_pred <- ms$enz_pred
  
  data <- mkdata_multi(data_path,0.1,0.01,10,grp,FALSE)
  num_g <- length(grp)
  enz <- data$enz
  num_sub <- ceiling(sqrt(num_rc))

  pdf(paste0(save_dir, "/hist_enz_pred.pdf"),width=10,height=10)
  par(mfrow=c(num_sub,num_sub)) 
  for(r in 1:num_rc){
    for(g in 1:num_g){
      enz_pred_now <- enz_pred[,r,data$idx_g[g,1]:data$idx_g[g,2]]
      enz_now <- enz[r,data$idx_g[g,1]:data$idx_g[g,2]]
      if(g==1){
        hist(enz_pred_now,
             col = paste0(col_list[g],"40"),
             border = col_list[g],
             main=rxn_names[r],
             xlim = c(0,2),
             xlab = "")
      }else{
        hist(enz_pred_now,
             col = paste0(col_list[g],"40"),
             border = col_list[g],
             main=rxn_names[r],
             xlim = c(0,2),
             axes = FALSE,
             add = TRUE)
      }
      rug(enz_now,ticksize=0.03,side=1,lwd=1,col=col_list[g])
    }
  }
  dev.off()
}

check_enz_fitting2 <- function(ms,grp,data_path,col_list,save_dir){
  enz_pred <- ms$enz_pred2
  data <- mkdata_multi(data_path,0.1,0.01,10,grp,FALSE)
  num_g <- length(grp)
  enz <- data$enz
  
  int_list_tmp <- mkdata_char(data_path,"/int_list.txt")
  num_rc <- as.integer(int_list_tmp[int_list_tmp[,1]=="num_rc",2])
  rxn_names <- mkdata_vec(data_path,"/rxn_names_include.txt")
  num_sub <- ceiling(sqrt(num_rc))
  
  pdf(paste0(save_dir, "/hist_enz_pred2.pdf"),width=10,height=10)
  par(mfrow=c(num_sub,num_sub)) 
  for(r in 1:num_rc){
    for(g in 1:num_g){
      enz_pred_now <- enz_pred[,r,data$idx_g[g,1]:data$idx_g[g,2]]
      enz_now <- enz[r,data$idx_g[g,1]:data$idx_g[g,2]]
      if(g==1){
        hist(enz_pred_now,
             col = paste0(col_list[g],"40"),
             border = col_list[g],
             main=rxn_names[r],
             xlim = c(0,2),
             xlab = "")
      }else{
        hist(enz_pred_now,
             col = paste0(col_list[g],"40"),
             border = col_list[g],
             main=rxn_names[r],
             xlim = c(0,2),
             axes = FALSE,
             add = TRUE)
      }
      rug(enz_now,ticksize=0.03,side=1,lwd=1,col=col_list[g])
    }
  }
  dev.off()
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

dens_flux_all <- function(ms,data_path,col_list1,col_list2,save_dir){
    ms_now_tmp <- ms["v"][[1]]
    dim_tmp <- dim(ms_now_tmp)
    ms_now_tmp2 <- ms_now_tmp[,,1]+ms_now_tmp[,,2]
    ms_now <- abind::abind(ms_now_tmp,ms_now_tmp2,along=3)
    min_now <- floor(min(ms_now))
    max_now <- ceiling(max(ms_now))
    idx_include <- mkdata_vec(data_path,"/idx_include.txt")

    rxn_names <- mkdata_vec(data_path,"/rxn_names_include.txt")
    rxn_names <- append(rxn_names,"G6pc")
    num_sub <- length(rxn_names)
    
    pdf(paste0(save_dir, "/dens_v_all.pdf"),width=5,height=num_sub*1.5)
    par(mar = c(2,4,2,2))#lower,left,upper,right
    par(mfrow=c(num_sub,1))
    # if parameters are 3 dimensions
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
    dev.off()
}

dens_flux_prior_all <- function(ms,data_path,col_list1,col_list2,save_dir){
  mu_vi <- ms["mu_vi"][[1]]
  iter <- dim(mu_vi)[1]
  # mu_v
  data <- mkdata_multi(data_path,0.1,0.01,10,grp,FALSE)
  idx_calc <- c(1,16,2,3,4,5,17,6,7,8,9,10,11,18,19,20,21,12,13,14,15,22)
  ms_now <- array(0,dim=c(iter,data$num_g,data$num_r))
  for(i in 1:iter){
    tmp <- t(data$Nu %*% t(mu_vi[i,,]))
    ms_now[i,,] <- tmp[,idx_calc]
  }
  
  min_now <- floor(min(ms_now))
  max_now <- ceiling(max(ms_now))
  rxn_names <- mkdata_vec(data_path,"/rxn_names.txt")
  num_sub <- length(rxn_names)
  
  pdf(paste0(save_dir, "/dens_mu_vi_all.pdf"),width=5,height=num_sub*1.5)
  par(mar = c(2,4,2,2))#lower,left,upper,right
  par(mfrow=c(num_sub,1))
  # if parameters are 3 dimensions
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
  dev.off()
}

trace_param <- function(fit,save_dir){
  ms <- rstan::extract(fit)
  param_list <- names(ms)
  for(par in param_list){
    if(par != "log_lik" && par != "y" && par !="x" && par != "y_eff" &&
       par != "enz_pred" && par != "enz_pred2"
       ){
      ### traceplot
      pdf(paste0(save_dir, "/traceplot_", par ,".pdf"),width=dim(ms)[2],height=dim(ms)[3])
      plot_now <- traceplot(fit, pars=par,inc_warmup=T)
      print(plot_now)
      dev.off()
    }
  }
}

dens_param <- function(pars_plot,ms,data_path,fname,save_dir){
  # num_sub <- length(pars_plot)
  num_sub <- 0
  for(par in pars_plot){
    ms_now <- ms[par][[1]]
    if(length(dim(ms_now))==1){
      num_sub <- num_sub+1
    }else{
      num_sub <- num_sub + ncol(ms_now)
    }
  }
  pdf(paste0(save_dir, "/dens_", fname, ".pdf"),width=5,height=num_sub*1.5)
  par(mar = c(2,4,2,2))#lower,left,upper,right
  par(mfrow=c(num_sub,1))
  for(par in pars_plot){
    ms_now <- ms[par][[1]]
    min_now <- floor(min(ms_now))
    max_now <- ceiling(max(ms_now))
    if(length(dim(ms_now))==1){
      d <- density(ms_now)
      m_now <- median(ms_now)
      sd_now <- sd(ms_now)
      main_now <- paste0(par," (median=", round(m_now,2), ", sd=", round(sd_now,2),")" )
      plot(d,xlim=c(min_now,max_now),main=main_now,xlab="",cex=5)
      polygon(d,col="grey")
    }else{
      for(i in 1:dim(ms_now)[2]){
        d <- density(ms_now[,i])
        m_now <- median(ms_now[,i])
        sd_now <- sd(ms_now[,i])
        main_now <- paste0(par," (median=", round(m_now,2), ", sd=", round(sd_now,2),")" )
        plot(d,xlim=c(min_now,max_now),main=paste0(main_now,i),xlab="",cex=5)
        polygon(d,col="grey")
      }
    }
    
  }
  dev.off()
}

output_summary <- function(fit,save_dir){
  tmp <- summary(fit)$summary
  write.csv(tmp,paste0(save_dir,'/summary_out.csv'))
}

flux_pca <- function(data_path,data,ms,save_dir){
  ######################
  v_pred <- ms$v
  
  num_g <- dim(v_pred)[2]
  v_data <- array(0,dim=c(20,num_g,50))
  v_data[,1,] <- mkdata(data_path,"/flux_smpl.txt")
  v_data[,2,] <- mkdata(data_path,"/flux_mut_smpl1.txt")
  v_data[,3,] <- mkdata(data_path,"/flux_mut_smpl2.txt")
  v_data[,4,] <- mkdata(data_path,"/flux_mut_smpl3.txt")
  v_data[,5,] <- mkdata(data_path,"/flux_mut_smpl4.txt")
  
  iter <- dim(v_pred)[1]
  v <- array(0,dim=c((iter+50)*num_g,20))
  col <- c("black","grey","red","pink","blue","lightblue","green","lightgreen")
  col <- c("#000000","#90909040","#ff0000","#ff00ff40","#0000ff","#00ffff40","#008000","#00ff0040","#ff8000","#ff800040")
  col_mat <- rep(0,(iter+50)*num_g)
  for(r in 1:20){
    for(g in 1:num_g){
      idx_s <- 50*(g-1)+1
      idx_e <- 50*g
      v[idx_s:idx_e,r] <- c(v_data[r,g,])
      col_mat[idx_s:idx_e] <- col[g*2-1]
      idx_s <- 50*4 + iter*(g-1) + 1
      idx_e <- 50*4 + iter*g
      v[idx_s:idx_e,r] <- c(v_pred[,g,r])
      col_mat[idx_s:idx_e] <- col[g*2]
    }
  }
  
  rpca <- prcomp(x=v,scale=T)
  
  sink(paste0(save_dir,"/summary_rpca.txt"))
  tmp <- summary(rpca)$importance
  print(tmp)
  sink()
  # PC1      PC2      PC3       PC4       PC5       PC6       PC7
  # Standard deviation     3.934516 1.706399 1.076633 0.5006535 0.3756614 0.1547765 0.1045994
  # Proportion of Variance 0.774020 0.145590 0.057960 0.0125300 0.0070600 0.0012000 0.0005500
  # Cumulative Proportion  0.774020 0.919610 0.977570 0.9901000 0.9971600 0.9983500 0.9989000
  
  max_pc1 <- ceiling(max(rpca$x[,1]))
  min_pc1 <- floor(min(rpca$x[,1]))
  max_pc2 <- ceiling(max(rpca$x[,2]))
  min_pc2 <- floor(min(rpca$x[,2]))
  max_pc3 <- ceiling(max(rpca$x[,3]))
  min_pc3 <- floor(min(rpca$x[,3]))
  
  pdf(paste0(save_dir, "/rpca_pc1_pc2.pdf"),width=8,height=8)
  par(mar = c(2,4,2,2))#lower,left,upper,right
  par(las=1)
  plot(x=NULL,type="n",xlab="PC1",ylab="PC2",xlim=c(min_pc1,max_pc1),ylim=c(min_pc2,max_pc2),xaxs="i",yaxs="i",xaxt="n",yaxt="n",bty="n")
  axis(side=1,at=seq(min_pc1,max_pc1,1),tck=1.0,lty="dotted",lwd=0.5,col="#dddddd")
  axis(side=2,at=seq(min_pc2,max_pc2,1),tck=1.0,lty="dotted",lwd=0.5,col="#dddddd")
  points(x=rpca$x[,1],y=rpca$x[,2],pch=20,col=col_mat)
  dev.off()
  
  pdf(paste0(save_dir, "/rpca_pc1_pc3.pdf"),width=8,height=8)
  par(mar = c(2,4,2,2))#lower,left,upper,right
  par(las=1)
  plot(x=NULL,type="n",xlab="PC1",ylab="PC3",xlim=c(min_pc1,max_pc1),ylim=c(min_pc3,max_pc3),xaxs="i",yaxs="i",xaxt="n",yaxt="n",bty="n")
  axis(side=1,at=seq(min_pc1,max_pc1,1),tck=1.0,lty="dotted",lwd=0.5,col="#dddddd")
  axis(side=2,at=seq(min_pc3,max_pc3,1),tck=1.0,lty="dotted",lwd=0.5,col="#dddddd")
  points(x=rpca$x[,1],y=rpca$x[,3],pch=16,col=col_mat)
  dev.off()
  
  library(scatterplot3d)
  pdf(paste0(save_dir, "/rpca_pc1_pc2_pc3.pdf"),width=8,height=8)
  par(mar = c(2,4,2,2))#lower,left,upper,right
  scatterplot3d(x=rpca$x[,1],y=rpca$x[,2],z=rpca$x[,3],
                color=col_mat)
  dev.off()
  
  save(data,file = paste0(save_dir,"/data.RData"))
  save(ms,file = paste0(save_dir,"/ms.RData"))
}
