functions{
  vector calc_mu(vector mu_vi,matrix Nu,int[] idx_calc,int num_r){
    vector[num_r] mu_v = Nu*(mu_vi);
    return mu_v[idx_calc];
  }
  
  matrix calc_Sigma2(vector mu_vi,matrix Nu,matrix Ne,matrix Sp,matrix Sm,real c_v,real c_e,int num_r,int num_m,int num_ri,int num_mc,int[] idx_calc){
    vector[num_ri] diag_mu_vi;
    vector[num_r] mu_v;
    vector[num_mc] dxdt;
    vector[num_mc] sigma_e;
    matrix[num_r,num_r] Sigma_tmp;
    for(r in 1:num_ri){
      diag_mu_vi[r] = (c_v*mu_vi[r])^2;
    }
    mu_v = Nu*mu_vi;
    dxdt = (fabs(Sp*mu_v)+fabs(Sm*mu_v))./2;
    for(m in 1:num_mc){
      sigma_e[m] = (c_e*dxdt[m])^2;
    }
    Sigma_tmp= Nu*(diag_matrix(diag_mu_vi))*Nu' + Ne*diag_matrix(sigma_e)*Ne';
    return Sigma_tmp[idx_calc,idx_calc];
  }
  
  real calc_log_prior(vector mu_vi,matrix Nu,matrix Ne,matrix Sp,matrix Sm,real c_v,real c_e,int num_r,int num_m,int num_ri,int num_mc,
                      int[] idx_calc,int num_rc, vector v){
   matrix[num_rc,num_rc] Sigma_v;
   vector[num_rc] mu_v;
   real target_prior;
   mu_v = calc_mu(mu_vi,Nu,idx_calc,num_r);
   Sigma_v = calc_Sigma2(mu_vi,Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc);
   target_prior = multi_normal_lpdf(v | mu_v, Sigma_v);
   return target_prior;
 }
  vector[] calc_x(vector mean_v,vector[] sub,vector[] pro,vector[] met_eff,
                real[] a,real[] b,real[] e_cofactor,real[] e_allo,
                vector c_OAA, vector c_Pyr,
                int num_rc,int num_s){
    vector[num_s] x[num_rc];
    vector[num_s] ones = rep_vector(1.0,num_s);

    // Pgm2
    x[1] = ones ./(mean_v[1]*(1+a[1]*sub[1]-b[1]*pro[1]));
    // Gpi1
    x[2] = ones ./(mean_v[2]*(1+a[2]*sub[2]-b[2]*pro[2]));
    // Fbp1 fbp1_i_amp fbp1_a_cit
    x[3] = ones ./(mean_v[3]*(1+a[3]*sub[3]-b[3]*pro[3]
          -e_allo[1]*met_eff[8] + e_allo[2]*met_eff[9]));
    // Gpd1, e_nadh[1:2]
    x[4] = ones ./(mean_v[4]*(1+a[4]*sub[4]-b[4]*pro[4]
          +e_cofactor[1]*met_eff[1]-e_cofactor[2]*met_eff[2]));
    // Pgam1
    x[5] = ones ./(mean_v[5]*(1+a[5]*sub[5]-b[5]*pro[5]));
    // Eno1
    x[6] = ones ./(mean_v[6]*(1+a[6]*sub[6]-b[6]*pro[6]));
    // Pklr e_atp[1:2] pklr_a_f16p pklr_i_ala pklr_i_phe
    x[7] = ones ./(mean_v[7]*(1+a[7]*sub[7]
          +e_cofactor[9]*met_eff[4]
          -e_cofactor[10]*met_eff[3]
          +e_allo[3]*met_eff[10]
          -e_allo[4]*met_eff[11]
          -e_allo[5]*met_eff[13]));
    // Pck1 e_gtp[1]
    x[8] = ones ./(mean_v[8]*(1+a[8]*c_OAA
          +e_cofactor[12]*met_eff[6]));
    //Ldha e_nadh[3:4]
    x[9] = ones ./(mean_v[9]*(1+a[9]*sub[9]-b[9]*c_Pyr
          +e_cofactor[3]*met_eff[1]-e_cofactor[4]*met_eff[2]));
    // Gpt
    x[10] = ones ./(mean_v[10]*(1+a[10]*sub[10]-b[10]*c_Pyr));
    // Pcx pcx_a_accoa e_atp[3]
    x[11] = ones ./(mean_v[11]*(1+a[11]*c_Pyr
            +e_cofactor[11]*met_eff[3]
            +e_allo[6]*met_eff[13]));
    //Cs
    x[12] = ones ./(mean_v[12]*(1+a[12]*sub[12]+b[12]*c_OAA));
    //Sdha
    x[13] = ones ./(mean_v[13]*(1+a[13]*sub[13]-b[13]*pro[13]
            +e_cofactor[13]*met_eff[7]));
    // Fh1
    x[14] = ones ./(mean_v[14]*(1+a[14]*sub[14]-b[14]*pro[14]));
    //Mdh2 e_nadh[5:6]
    x[15] = ones ./(mean_v[15]*(1+a[15]*sub[15]-b[15]*c_OAA
            +e_cofactor[5]*met_eff[1]-e_cofactor[6]*met_eff[2]));
    //Glud1 e_nadh[7:8] glud1_i_atp glud1_a_adp glud1_i_gtp glud1_i_leu
    x[16] = ones ./(mean_v[16]*(1+a[16]*sub[16]
            +e_cofactor[7]*met_eff[1]-e_cofactor[8]*met_eff[2]
            -e_allo[7]*met_eff[3]+e_allo[8]*met_eff[4]
            -e_allo[9]*met_eff[6]-e_allo[10]*met_eff[14]));
    
    return x;
 }
  vector prep_mu_vi_ref(vector mu_vi_wt_tmp,real v_Pcx,real v_Cs,int num_ri){
    vector[num_ri] mu_vi_ref;
    mu_vi_ref[1:4] = mu_vi_wt_tmp[1:4];
    mu_vi_ref[5] = v_Pcx;
    mu_vi_ref[6] = v_Cs;
    mu_vi_ref[7] = 1-sum(mu_vi_wt_tmp[1:4]);
    return mu_vi_ref;
  } 
  
  vector prep_mean_v(vector mu_vi_wt_tmp,real v_Pcx,real v_Cs,vector[] mu_vi_g,
                     int num_r,int num_ri,int num_rc,int num_g,
                     matrix Nu,int[] idx_calc,int[,] idx_g){
    vector[num_rc] mean_v;
    vector[num_ri] mu_vi_ref;
    vector[num_rc] mu_v[num_g];
    real sum_v[num_rc] = rep_array(0.0,num_rc);
    
    // mu_vi for reference group(WT)
    mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,v_Pcx,v_Cs,num_ri);
    mu_v[1] = calc_mu(mu_vi_ref,Nu,idx_calc,num_r);
    for(g in 2:num_g){
      mu_v[g] = calc_mu(mu_vi_g[g-1],Nu,idx_calc,num_r);
    }
    for(g in 1:num_g){
      for(r in 1:num_rc){
        sum_v[r] += mu_v[g][r]*(idx_g[g,2]-idx_g[g,1]+1);
      }
    }
    for(r in 1:num_rc){
        mean_v[r] = sum_v[r]/idx_g[num_g,2];
    }
    return to_vector(mean_v);
  }
  
  vector[] prep_r_p_all(vector[] r_p,vector r_p_tmp,
                      int num_rc,int num_g,int num_p,
                      int[] idx_p,int[] idx_p_,int[,] idx_g,
                      vector[] rna){
    vector[num_rc] r_p_all[num_g];
    real sum_r[num_p] = rep_array(0.0,num_p);
    
    for(g in 1:(num_g-1)){
      // prepare r_p_all from r_p and r_p_tmp
      r_p_all[g] = r_p[g];
      for(p in 1:num_p){//Gpt,Glud1
        sum_r[p] += r_p[g][idx_p[p]]*sum(rna[idx_p[p]][idx_g[g,1]:idx_g[g,2]]);
      }
    }
    // prepare r_p_all from r_p and r_p_tmp
    r_p_all[num_g][idx_p_] = r_p_tmp;
    for(p in 1:num_p){
      r_p_all[num_g][idx_p[p]] = -sum_r[p] / sum( rna[idx_p[p]][idx_g[num_g,1]:idx_g[num_g,2]] );
    }
    return r_p_all;
                      }
  vector[] prep_y(vector[] r_p,vector r_p_tmp, vector[] rna, vector[] v_x,
                  int[,] idx_g,int[] idx_p,int[] idx_p_,int num_g,int num_rc,int num_p){
    vector[idx_g[num_g,2]] y[num_rc];
    vector[num_rc] r_p_all[num_g];
    
    r_p_all = prep_r_p_all(r_p,r_p_tmp,
                             num_rc,num_g,num_p,
                             idx_p,idx_p_,idx_g,
                             rna);
    for(g in 1:num_g){
      for(r in 1:num_rc){
        y[r][idx_g[g,1]:idx_g[g,2]] = v_x[r][idx_g[g,1]:idx_g[g,2]] ./ (1+r_p_all[g][r]);
      }
    }
    return y;
  }
  
  vector[] prep_OAA_Pyr(real[] c_OAA,real[] c_Pyr,
                        int[,] idx_g,int num_g){
    vector[idx_g[num_g,2]] c_out[2];
    real sum_c_OAA = 0.0;
    real sum_c_Pyr = 0.0;
    int num_smpl;
    
    num_smpl = idx_g[num_g,2];
    for(g in 1:(num_g-1)){
      // sum used to adjust c_OAA and c_Pyr to make sure the mean equals 1
      sum_c_OAA += c_OAA[g]*(idx_g[g,2]-idx_g[g,1]+1);
      sum_c_Pyr += c_Pyr[g]*(idx_g[g,2]-idx_g[g,1]+1);
      // make vector for c_OAA,c_Pyr
      c_out[1][idx_g[g,1]:idx_g[g,2]] = rep_vector(log(c_OAA[g]),idx_g[g,2]-idx_g[g,1]+1);
      c_out[2][idx_g[g,1]:idx_g[g,2]] = rep_vector(log(c_Pyr[g]),idx_g[g,2]-idx_g[g,1]+1);
    }
    // make vector for c_OAA,c_Pyr
    c_out[1][idx_g[num_g,1]:idx_g[num_g,2]] = rep_vector(log( (num_smpl-sum_c_OAA) / (idx_g[num_g,2]-idx_g[num_g,1]+1) ),
                                                        idx_g[num_g,2]-idx_g[num_g,1]+1);
    c_out[2][idx_g[num_g,1]:idx_g[num_g,2]] = rep_vector(log( (num_smpl-sum_c_Pyr) / (idx_g[num_g,2]-idx_g[num_g,1]+1) ),
                                                        idx_g[num_g,2]-idx_g[num_g,1]+1);
    return c_out;
  }
}

data{
    int<lower=1> num_r;// # reaction in the pathway
    int<lower=1> num_ri;// # independent flux
    int<lower=1> num_rd;// # dependent flux
    int<lower=1> num_rc;// # reaction flux for kinetic modeling
    int<lower=1> num_m;// # metabolites
    int<lower=1> num_mc;// # not end metabolites
    int<lower=1> num_b;// # number of b
    int<lower=1> num_ri_wt;// # number of independent flux in reference(WT)
    int<lower=1> num_smpl;// # number of sample (including all groups)
    int<lower=1> num_g;// # group
    int<lower=1> num_p;// # reactions for enzyme prediction
    matrix[num_m,num_r] S;// stoichiometry matrix
    matrix[num_r,num_ri] N;// null space
    matrix[num_r,num_ri] Nu;
    matrix[num_r,num_rd] Ne;
    matrix[num_mc,num_r] Sp;
    matrix[num_mc,num_r] Sm;
    int idx_calc[num_rc];// indec for kinetic modeling
    vector[num_smpl] enz[num_rc];// enzyme data
    vector[num_smpl] sub[num_rc];// substrate data
    vector[num_smpl] pro[num_rc];// product data
    vector[num_smpl] rna[num_rc];// transcript data
    int idx_g[num_g,2];// index(start,end) list for group
    int idx_p[num_p];// index for enzyme prediction
    int idx_p_[num_rc-num_p];//index for NOT enzyme prediction
    int idx_b[num_b-1];// index for b
    int is_irrev[num_rc];// logical vector 0(reversible) or 1(irreversible)
    real<lower=0,upper=1> c_v;// coefficient of variance for v
    real<lower=0,upper=1> c_e;// coefficient of variance for \dot{x}
    vector[num_smpl] enz_eff[1];// enzymes consisting complex, Sdhb
    vector[num_smpl] rna_eff[1];// transcripts of enzymes consisting complex, Sdhb
    
    vector[num_smpl] met_eff[14];// metabolite data for effectors(cofactors and allosteric effectors)
    // 1"NAD+" 2"NADH" 3"ATP" 4"ADP" 5"GDP"" 6"GTP" 7"FAD" 8"AMP" 9"Cit"
    // 10"F1,6P" 11"Ala" 12"AcCoA" 13"Phe" 14"Leu"
    
    int v_max;// flux range
}
parameters{
    vector[num_rc] v[num_g];// flux
    
    real<lower=0> a[num_rc];// elasticity coefficient for substrate
    real<lower=0> b[num_b-1];// elasticity coefficient for products
    
    // cofactors
    real<lower=0> e_cofactor[13];
    // real<lower=0,upper=1> e_nadh[8];
    // real<lower=0,upper=1> e_atp[3];
    // real<lower=0,upper=1> e_gtp;
    // real<lower=0,upper=1> e_fad;
    
    // allosteric effectors
    real<lower=0> e_allo[10];
    // 1 real<lower=0,upper=1> fbp1_i_amp;//met_eff[8]
    // 2 real<lower=0,upper=1> fbp1_a_cit;//met_eff[9]
    // 3 real<lower=0,upper=1> pklr_a_f16p;//met_eff[10]
    // 4 real<lower=0,upper=1> pklr_i_ala;//met_eff[12]
    // 5 real<lower=0,upper=1> pklr_i_phe;
    // 6 real<lower=0,upper=1> pcx_a_accoa;//met_eff[13]
    // 7 real<lower=0,upper=1> glud1_i_atp;
    // 8 real<lower=0,upper=1> glud1_a_adp;
    // 9 real<lower=0,upper=1> glud1_i_gtp;
    // 10 real<lower=0,upper=1> glud1_i_leu;
    
    real<lower=0,upper=4> c_OAA[num_g-1];// OAA concentration relative to reference state(WT)
    real<lower=0,upper=4> c_Pyr[num_g-1];// pyruvate concentration relative to reference state(WT)
    
    // mu_vi for reference group
    vector<lower=-v_max,upper=v_max>[4] mu_vi_wt_tmp;
    real<lower=-v_max,upper=v_max> v_Pcx;
    real<lower=-v_max,upper=v_max> v_Cs;
    
    vector<lower=-v_max,upper=v_max>[num_ri] mu_vi_g[num_g-1];// mu_vi for other group
    
    real<lower=0> sigma_n;// error term for likelihood of metabolic submodel
    real<lower=0> sigma_n2;// error term for likelihood of protein expression submodel
    
    vector[num_rc] r_p[num_g-1];// protein turnover coefficient (except adjusted group)
    vector[num_rc-num_p] r_p_tmp;// protein turnover coefficient (adjusted group)
    
    real r_p_eff[num_g];// protein turnover coefficient for Sdhb
    real<lower=0> sigma_p;// error term for prior in protein expression submodel
    
}
transformed parameters{
  vector[num_smpl] x[num_rc];
  vector<lower=0>[num_smpl] y_eff;
  vector<lower=0>[num_smpl] y[num_rc];
  {
      real b2[num_rc]  = rep_array(0.0,num_rc);
      vector[num_rc] mean_v;
      vector[num_smpl] c_out[2];
      // real e_cofactor[10]=rep_array(0.0,10);
      // real e_allo[10]=rep_array(0.0,10);
      vector[num_smpl] v_x[num_rc];

      for(g in 1:num_g){
        // prepare y and enz
        y_eff[idx_g[g,1]:idx_g[g,2]] = enz_eff[1][idx_g[g,1]:idx_g[g,2]] ./ (1+r_p_eff[g]);
      }
      
      // calculate mean of all the prior mu_v
      mean_v = prep_mean_v(mu_vi_wt_tmp,v_Pcx,v_Cs,mu_vi_g,
                           num_r,num_ri,num_rc,num_g,
                           Nu,idx_calc,idx_g);
      b2[idx_b] = b;
      
      c_out = prep_OAA_Pyr(c_OAA,c_Pyr,idx_g,num_g);
      // calculate X
      x = calc_x(mean_v,sub,pro,met_eff,
                 a,b2,e_cofactor,e_allo,
                 c_out[1],c_out[2],
                 num_rc,num_smpl);
                 
      for(r in 1:num_rc){
        for(g in 1:num_g){
          v_x[r][idx_g[g,1]:idx_g[g,2]] = v[g][r] * x[r][idx_g[g,1]:idx_g[g,2]];
        }
      }
      
      y = prep_y(r_p,r_p_tmp,rna,v_x,
                 idx_g,idx_p,idx_p_,num_g,num_rc,num_p);
    }
}
model{
  vector[num_smpl] v_x[num_rc];
  vector[num_ri] mu_vi_ref;
  
  for(r in 1:num_rc){
    for(g in 1:num_g){
      v_x[r][idx_g[g,1]:idx_g[g,2]] = v[g][r] * x[r][idx_g[g,1]:idx_g[g,2]];
    }
  }
  
  for(r in 1:num_rc){
      if(r==idx_p[1]){
        target += normal_lpdf(rna[r] | y[r], sqrt(sigma_n2));
      }else if(r==idx_p[2]){
        target += normal_lpdf(rna[r] | y[r], sqrt(sigma_n2));
      }else if(r==13){
        target += normal_lpdf(enz[13] .* enz_eff[1] | v_x[r], sqrt(sigma_n) );
        target += normal_lpdf(rna[r] | y[r], sqrt(sigma_n2));
        target += normal_lpdf(rna_eff[1] | y_eff, sqrt(sigma_n2));
      }else{
        target += normal_lpdf(enz[r] | v_x[r], sqrt(sigma_n) );
        target += normal_lpdf(rna[r] | y[r], sqrt(sigma_n2));
      }
    }
    
  mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,v_Pcx,v_Cs,num_ri);
  target += calc_log_prior(mu_vi_ref,Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[1]);
  for(g in 2:num_g){
    target += calc_log_prior(mu_vi_g[g-1],Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[g]);
  }
  
  for(g in 1:num_g){
    if(g != num_g){
      target += normal_lpdf(r_p[g]|0,sigma_p);
    }else{
      target += normal_lpdf(r_p_tmp[1]|0,sigma_p);
    }
    target += normal_lpdf(r_p_eff[g]|0,sigma_p);
  }

  target += normal_lpdf(a|0,1);
  target += normal_lpdf(b|0,1);
  target += normal_lpdf(e_cofactor|0,1);
  target += normal_lpdf(e_allo|0,1);
  
  target += cauchy_lpdf(sigma_p|0,0.5);
  target += cauchy_lpdf(sigma_n|0,0.5);
  target += cauchy_lpdf(sigma_n2|0,0.5);
    
}
generated quantities{
  vector[num_smpl] log_lik[num_rc];
  real log_prior[num_g];
  vector[num_smpl] enz_pred[num_rc];
  vector[num_smpl] rna_pred[num_rc];
  {
    vector[num_smpl] enz_Sdha;
    vector[num_ri] mu_vi_ref;
    
    enz_Sdha = enz[13] .* enz_eff[1];
    mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,v_Pcx,v_Cs,num_ri);
    
    log_prior[1] = calc_log_prior(mu_vi_ref,Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[1]);
    for(g in 2:num_g){
      log_prior[g] = calc_log_prior(mu_vi_g[g-1],Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[g]);
    }
    
    for(g in 1:num_g){
      for(r in 1:num_rc){
        for(n in idx_g[g,1]:idx_g[g,2]){
          if(r==idx_p[1]){
            log_lik[r][n] = normal_lpdf(rna[r][n] | y[r][n], sqrt(sigma_n2) );
          }else if(r==idx_p[2]){
            log_lik[r][n] = normal_lpdf(rna[r][n] | y[r][n], sqrt(sigma_n2) );
          }else if(r==13){
            log_lik[r][n] = normal_lpdf(enz_Sdha[n] | v[g][r]*x[r][n], sqrt(sigma_n) )
                        + normal_lpdf(rna[13][n] | y[r][n], sqrt(sigma_n2) )
                        + normal_lpdf(rna_eff[1][n] | y_eff[n], sqrt(sigma_n2));
          }else{
            log_lik[r][n] = normal_lpdf(enz[r][n] | v[g][r]*x[r][n], sqrt(sigma_n) )
                        + normal_lpdf(rna[r][n] | y[r][n], sqrt(sigma_n2) );
          }
          enz_pred[r][n] = normal_rng(v[g][r]*x[r][n], sqrt(sigma_n));
          rna_pred[r][n] = normal_rng(y[r][n], sqrt(sigma_n2));
        }
      }
    }
  }
}
