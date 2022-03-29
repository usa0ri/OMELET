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
    Sigma_tmp= Nu*(diag_matrix(diag_mu_vi))*Nu'+ Ne*diag_matrix(sigma_e)*Ne';
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

  vector[] calc_x(vector mean_v,vector[] sub,vector[] pro,
                real[] a,real[] b,
                vector[] met_eff,real[] e_effs,
                vector[] met_est,
                int num_rc,int num_s){
    vector[num_s] x[num_rc];
    vector[num_s] ones = rep_vector(1.0,num_s);

    // Pgm2
    x[1] = ones ./(mean_v[1]*(1+a[1]*sub[1]-b[1]*pro[1]));
    // Gpi1
    x[2] = ones ./(mean_v[2]*(1+a[2]*sub[2]-b[2]*pro[2]));
    // Fbp1
    x[3] = ones ./(mean_v[3]*(1+a[3]*sub[3]-b[3]*pro[3]
          -e_effs[1]*met_eff[8]+e_effs[2]*met_eff[9]));
    // Gpd1
    x[4] = ones ./(mean_v[4]*(1+a[4]*sub[4]-b[4]*pro[4]
          +e_effs[3]*met_eff[1]-e_effs[4]*met_eff[2]));
    // Pgam1
    x[5] = ones ./(mean_v[5]*(1+a[5]*sub[5]-b[5]*pro[5]));
    // Eno1
    x[6] = ones ./(mean_v[6]*(1+a[6]*sub[6]-b[6]*pro[6]));
    // Pklr
    x[7] = ones ./(mean_v[7]*(1+a[7]*sub[7]-b[7]*pro[7]
          +e_effs[5]*met_eff[4]+e_effs[6]*met_eff[10]-e_effs[7]*met_eff[11]-e_effs[8]*met_eff[13]));
    // Pck1
    x[8] = ones ./(mean_v[8]*(1+a[8]*met_est[1]-b[8]*pro[8]
          +e_effs[9]*met_eff[6]));
    // Ldha
    x[9] = ones ./(mean_v[9]*(1+a[9]*sub[9]-b[9]*met_est[2]
          +e_effs[10]*met_eff[1]-e_effs[11]*met_eff[2]));
    // Gpt
    x[10] = ones ./(mean_v[10]*(1+a[10]*sub[10]-b[10]*met_est[2]));
    // Pcx
    x[11] = ones ./(mean_v[11]*(1+a[11]*met_est[2]-b[11]*pro[11]
          +e_effs[12]*met_eff[3]+e_effs[13]*met_eff[12]));
    // Cs
    x[12] = ones ./(mean_v[12]*(1+a[12]*met_est[1]-b[12]*pro[12]));
    // Sdha
    x[13] = ones ./(mean_v[13]*(1+a[13]*sub[13]-b[13]*pro[13]
          +e_effs[14]*met_eff[7]));
    // Fh1
    x[14] = ones ./(mean_v[14]*(1+a[14]*sub[14]-b[14]*pro[14]));
    // Mdh2
    x[15] = ones ./(mean_v[15]*(1+a[15]*sub[15]-b[15]*met_est[1]
          +e_effs[15]*met_eff[1]-e_effs[16]*met_eff[2]));
    // Glud1
    x[16] = ones ./(mean_v[16]*(1+a[16]*sub[16]-b[16]*pro[16]
          +e_effs[17]*met_eff[1]-e_effs[18]*met_eff[2]-e_effs[19]*met_eff[3]+e_effs[20]*met_eff[4]-e_effs[21]*met_eff[6]-e_effs[22]*met_eff[14]));
    return x;
 }
  vector prep_mu_vi_ref(vector mu_vi_wt_tmp,int num_ri, matrix Nu){
    vector[num_ri] mu_vi_ref = rep_vector(0,num_ri);
    mu_vi_ref[1:6] = mu_vi_wt_tmp[1:6];
    mu_vi_ref[7] = 1-Nu[1,]*mu_vi_ref;
    return mu_vi_ref;
 }
  vector prep_mean_v(vector mu_vi_wt_tmp,vector[] mu_vi_g,
                     int num_r,int num_ri,int num_rc,int num_g,
                     matrix Nu,int[] idx_calc,int[,] idx_g){
    vector[num_rc] mean_v;
    vector[num_ri] mu_vi_ref;
    vector[num_rc] mu_v[num_g];
    real sum_v[num_rc] = rep_array(0.0,num_rc);

    // mu_vi for reference group(WT)
    mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,num_ri,Nu);
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

  vector[] prep_r_p_all(vector[] r_p,vector r_p_tmp,int num_rc,int num_g,int num_p,int[] idx_p,int[] idx_p_,int[,] idx_g,vector[] rna){
    vector[num_rc] r_p_all[num_g];
    real sum_r[num_p] = rep_array(0.0,num_p);

    for(g in 1:(num_g-1)){
      r_p_all[g] = r_p[g];
      for(p in 1:num_p){
        sum_r[p] += r_p[g][idx_p[p]]*sum(rna[idx_p[p]][idx_g[g,1]:idx_g[g,2]]);
      }
    }

    r_p_all[num_g][idx_p_] = r_p_tmp;
    for(p in 1:num_p){
      r_p_all[num_g][idx_p[p]] = -sum_r[p] / sum( rna[idx_p[p]][idx_g[num_g,1]:idx_g[num_g,2]] );
    }
    return r_p_all;
  }

  vector[] prep_y(vector[] r_p,vector r_p_tmp, vector[] rna,vector[] v_x,
                  int[,] idx_g,int[] idx_p,int[] idx_p_,int num_g,int num_rc,int num_p){
    vector[idx_g[num_g,2]] y[num_rc];
    vector[num_rc] r_p_all[num_g];

    r_p_all = prep_r_p_all(r_p,r_p_tmp,num_rc,num_g,num_p,idx_p,idx_p_,idx_g,rna);
    for(g in 1:num_g){
      for(r in 1:num_rc){
        y[r][idx_g[g,1]:idx_g[g,2]] = v_x[r][idx_g[g,1]:idx_g[g,2]] ./ (1+r_p_all[g][r]);
      }
    }
    return y;
  }

  vector[] prep_met_est(vector[] c_input,int[,] idx_g,int num_g,int num_smpl,int num_met_est){
    vector[num_smpl] c_out[num_met_est];
    real sum_c[num_met_est] = rep_array(0.0,num_met_est);

    for(g in 1:(num_g-1)){
      for(c in 1:num_met_est){
        sum_c[c] += c_input[c][g]*(idx_g[g,2]-idx_g[g,1]+1);
        c_out[c][idx_g[g,1]:idx_g[g,2]] = rep_vector(log(c_input[c][g]),idx_g[g,2]-idx_g[g,1]+1);
      }
    }
    for(c in 1:num_met_est){
      c_out[c][idx_g[num_g,1]:idx_g[num_g,2]] = rep_vector(log((num_smpl-sum_c[c])/(idx_g[num_g,2]-idx_g[num_g,1]+1)),idx_g[num_g,2]-idx_g[num_g,1]+1);
    }
    return c_out;
  }
}

data{
    int<lower=0> num_r;
    int<lower=0> num_ri;
    int<lower=0> num_rd;
    int<lower=0> num_rc;
    int<lower=0> num_m;
    int<lower=0> num_mc;
    int<lower=0> num_ri_wt;
    int<lower=0> num_b;
    int<lower=0> num_smpl;
    int<lower=0> num_g;
    int<lower=0> num_p;
    int<lower=0> num_met_eff;
    int<lower=0> num_met_eff_pairs;
    int<lower=0> num_enz_cmplx;
    int<lower=0> num_met_est;
    matrix[num_m,num_r] S;
    matrix[num_r,num_ri] N;
    matrix[num_r,num_ri] Nu;
    matrix[num_r,num_rd] Ne;
    matrix[num_rd,num_r] Sp;
    matrix[num_rd,num_r] Sm;
    int idx_calc[num_rc];
    int idx_g[num_g,2];
    int idx_p[num_p];
    int idx_p_[num_rc-num_p];
    int idx_b[num_b];
    int is_indflux[num_r];
    int idx_indflux[num_ri];
    int idx_include[num_rc];
    int is_irrev[num_rc];
    int idx_irrev[num_rc-num_b];
    int is_cmplx[num_rc];
    int<lower=0> idx_fixed;
    vector[num_smpl] enz[num_rc];
    vector[num_smpl] rna[num_rc];
    vector[num_smpl] sub[num_rc];
    vector[num_smpl] pro[num_rc];
    vector[num_smpl] met_eff[num_rc];
    real<lower=0,upper=1> c_v;
    real<lower=0,upper=1> c_e;
    int v_max;
}

parameters{
    vector[num_rc] v[num_g];// flux
    real<lower=0> a[num_rc];// elasticity coefficient for substrate
    real<lower=0> b[num_b];// elasticity coefficient for products
    vector<lower=-v_max,upper=v_max>[6] mu_vi_wt_tmp;
    vector<lower=-v_max,upper=v_max>[num_ri] mu_vi_g[num_g-1];// mu_vi for other group
    real<lower=0> sigma_n;// error term for likelihood of enzymes

    real<lower=0> sigma_p;
    real<lower=0> sigma_n2;
    vector[num_rc] r_p[num_g-1];
    vector[num_rc-num_p] r_p_tmp;
    real<lower=0> e_effs[num_met_eff_pairs];// elasticity coefficient for effectors (cofactors and allosteric effectors)
    vector<lower=0,upper=4>[num_g-1] met_est[num_met_est];
}

transformed parameters{
  vector[num_smpl] x[num_rc];
  vector<lower=0>[num_smpl] y[num_rc];
  {
      real b2[num_rc]  = rep_array(0.0,num_rc);
      vector[num_rc] mean_v;
      vector[num_smpl] c_out[num_met_est];
      vector[num_smpl] v_x[num_rc];

      mean_v = prep_mean_v(mu_vi_wt_tmp,mu_vi_g,num_r,num_ri,num_rc,num_g,Nu,idx_calc,idx_g);
      b2[idx_b] = b;
      c_out = prep_met_est(met_est,idx_g,num_g,num_smpl,num_met_est);
      x = calc_x(mean_v,sub,pro,a,b2,met_eff,e_effs,c_out,num_rc,num_smpl);

      for(r in 1:num_rc){
        for(g in 1:num_g){
          v_x[r][idx_g[g,1]:idx_g[g,2]] = v[g][r] * x[r][idx_g[g,1]:idx_g[g,2]];
        }
      }
     y = prep_y(r_p,r_p_tmp,rna,v_x,idx_g,idx_p,idx_p_,num_g,num_rc,num_p);
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
        target += normal_lpdf(rna[r] | y[r], sqrt(sigma_n2) );
      }else if(r==idx_p[2]){
        target += normal_lpdf(rna[r] | y[r], sqrt(sigma_n2) );
      }else {
        target += normal_lpdf(enz[r] | v_x[r], sqrt(sigma_n) );
        target += normal_lpdf(rna[r] | y[r], sqrt(sigma_n2));
      }
  }

  for(g in 1:num_g){
    if(g != num_g){
      target += normal_lpdf(r_p[g]|0,sigma_p);
    }else{
      for(i in 1:num_p){
        target += normal_lpdf(r_p_tmp[i]|0,sigma_p);
      }
    }
  }

  mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,num_ri,Nu);
   target += calc_log_prior(mu_vi_ref,Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[1]);
  for(g in 2:num_g){
    target += calc_log_prior(mu_vi_g[g-1],Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[g]);
  }

  target += normal_lpdf(a|0,1);
  target += normal_lpdf(b|0,1);
  target += normal_lpdf(e_effs|0,1);
  target += cauchy_lpdf(sigma_n|0,0.5);
  target += cauchy_lpdf(sigma_p|0,0.5);
  target += cauchy_lpdf(sigma_n2|0,0.5);
  
}
generated quantities{
  vector[num_smpl] log_lik[num_rc];
  real log_prior[num_g];
  vector[num_smpl] enz_pred[num_rc];
  vector[num_smpl] rna_pred[num_rc];
{
    vector[num_ri] mu_vi_ref;

    mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,num_ri,Nu);
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
             }else {
             log_lik[r][n] = normal_lpdf(enz[r] | v[g][r]*x[r][n], sqrt(sigma_n) )
                           + normal_lpdf(rna[r] | y[r], sqrt(sigma_n2));
             }

          enz_pred[r][n] = normal_rng(v[g][r]*x[r][n], sqrt(sigma_n));
          rna_pred[r][n] = normal_rng(y[r][n], sqrt(sigma_n2));
    }
   }
  }
  }
}
