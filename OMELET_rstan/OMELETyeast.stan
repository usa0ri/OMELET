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
  
  vector[] calc_x(vector mean_v,vector[] sub1,vector[] sub2,vector[] pro1,vector[] pro2,vector[] met_eff,
                 real[] a,real a_2,real[] b,real b_2,real[] e_atp,real[] e_nad,real[] e_allo,
                 int num_rc,int num_s){
    vector[num_s] x[num_rc];
    vector[num_s] ones = rep_vector(1.0,num_s);
    
      //GND NADPH/NADP e_nadp[1:2]
      // x[1] = ones ./(mean_v[1]*(1+a[1]*sub1[1]-b[1]*pro1[1]
      //           +e_nadp[1]*met_eff[2]-e_nadp[2]*met_eff[1]));
      // //RKI
      // x[2] = ones ./(mean_v[2]*(1+a[2]*sub1[2]-b[2]*pro1[2]));
      // //RPE
      // x[3] = ones ./(mean_v[3]*(1+a[3]*sub1[3]-b[3]*pro1[3]));
      // //SOL
      // x[4] = ones ./(mean_v[4]*(1+a[4]*sub1[4]-b[4]*pro1[4]));
      // //TAL a_2[1] b_2[1]
      // x[5] = ones ./(mean_v[5]*(1+a[5]*sub1[5]+a_2[1]*sub2[5]-b[5]*pro1[5]-b_2[1]*pro2[5]));
      // //TKL_E4PF6P a_2[2] b_2[2]
      // x[6] = ones ./(mean_v[6]*(1+a[6]*sub1[6]+a_2[2]*sub2[6]-b[6]*pro1[6]-b_2[2]*pro2[6]));
      // //TKL_R5PS7P a_2[3] b_2[3]
      // x[7] = ones ./(mean_v[7]*(1+a[7]*sub1[7]+a_2[3]*sub2[7]-b[7]*pro1[7]-b_2[3]*pro2[7]));
      // //ZWF NADP/NADPH e_nadp[3:4]
      // x[8] = ones ./(mean_v[8]*(1+a[8]*sub1[8]-b[8]*pro1[8]
      //           +e_nadp[3]*met_eff[1]-e_nadp[4]*met_eff[2]));
      // 
      //HXT
      x[1] = ones ./(mean_v[1]*(1+a[1]*sub1[1]-b[1]*pro1[1]));
      //HXK ATP/ADP e_atp[1:2] HXK_T6P(i) e_allo[1]
      x[2] = ones ./(mean_v[2]*(1+a[2]*sub1[2]-b[2]*pro1[2]
                +e_atp[1]*met_eff[4]-e_atp[2]*met_eff[3]
                -e_allo[1]*met_eff[7]));
      // TPS UDP a_2[1]
      x[3] = ones ./(mean_v[3]*(1+a[3]*sub1[3]+a_2*sub2[3]-b[3]*pro1[3]));
      // TPP
      x[4] = ones ./(mean_v[4]*(1+a[4]*sub1[4]));
      x[5] = ones ./(mean_v[5]*(1+a[5]*sub1[5]-b[5]*pro1[5]));
      x[6] = ones ./(mean_v[6]*(1+a[6]*sub1[6]-b[6]*pro1[6]));
      x[7] = ones ./(mean_v[7]*(1+a[7]*sub1[7]-b[7]*pro1[7]));
      // PFK ATP/ADP e_atp[3:4] PFK_F26bP(a),PFK_AMP(a) e_allo[2:3]
      x[8] = ones ./(mean_v[8]*(1+a[8]*sub1[8]-b[8]*pro1[8]
                +e_atp[3]*met_eff[4]-e_atp[4]*met_eff[3]
                +e_allo[2]+e_allo[3]*met_eff[9]));
      // FBA b_2
      x[9] = ones ./(mean_v[9]*(1+a[9]*sub1[9]-b[9]*pro1[9]+b_2*pro2[9]));
      x[10] = ones ./(mean_v[10]*(1+a[10]*sub1[10]-b[10]*pro1[10]));
      // GPD NAD/NADH e_nad[1:2] GPD_F16bP(i??) e_allo[4]
      x[11] = ones ./(mean_v[11]*(1+a[11]*sub1[11]-b[11]*pro1[11]
                 +e_nad[1]*met_eff[5]-e_nad[2]*met_eff[6]
                 -e_allo[4]*met_eff[10]));
      // GPP        
      x[12] = ones ./(mean_v[12]*(1+a[12]*sub1[12]));
      // TDH NAD/NADH e_nad[3:4]
      x[13] = ones ./(mean_v[13]*(1+a[13]*sub1[13]-b[13]*pro1[13]
                   +e_nad[3]*met_eff[5]-e_nad[4]*met_eff[6]));
      // PGK ADP/ATP e_atp[5:6]            
      x[14] = ones ./(mean_v[14]*(1+a[14]*sub1[14]-b[14]*pro1[14]
                   +e_atp[5]*met_eff[3]-e_atp[6]*met_eff[4]));
      x[15] = ones ./(mean_v[15]*(1+a[15]*sub1[15]-b[15]*pro1[15]));
      x[16] = ones ./(mean_v[16]*(1+a[16]*sub1[16]-b[16]*pro1[16]));
      // PYK ADP/ATP e_atp[7:8], PYK_F16bP(a) e_allo[5]
      x[17] = ones ./(mean_v[17]*(1+a[17]*sub1[17]-b[17]*pro1[17]
                  +e_atp[7]*met_eff[3]-e_atp[8]*met_eff[4]
                  +e_allo[5]*met_eff[10]));
      x[18] = ones ./(mean_v[18]*(1+a[18]*sub1[18]-b[18]*pro1[18]));
      // acetate_branch NAD/NADH e_nad[5:6]
      x[19] = ones ./(mean_v[19]*(1+a[19]*sub1[19]-b[19]*pro1[19]
      +e_nad[5]*met_eff[5]-e_nad[6]*met_eff[6]));
      // ADH NADH/NAD e_nad[7:8]
      x[20] = ones ./(mean_v[20]*(1+a[20]*sub1[20]-b[20]*pro1[20]
      +e_nad[7]*met_eff[6]-e_nad[8]*met_eff[5]));
    return x;
  }
  
  vector prep_mu_vi_ref(vector mu_vi_wt_tmp,int num_ri){
    vector[num_ri] mu_vi_ref;
    //1 HXT=1
      mu_vi_ref[2:num_ri] = mu_vi_wt_tmp;
      mu_vi_ref[1] = (1-0.5*sum(mu_vi_ref[2:4]))/2;
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
    mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,num_ri);
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
                      
  vector[] prep_y(vector[] r_p,vector[] v_x,
                  int[,] idx_g,int num_g,int num_rc){
    vector[idx_g[num_g,2]] y[num_rc];
    vector[num_rc] r_p_all[num_g];
    
    for(g in 1:num_g){
      for(r in 1:num_rc){
        y[r][idx_g[g,1]:idx_g[g,2]] = v_x[r][idx_g[g,1]:idx_g[g,2]] ./ (1+r_p[g][r]);
      }
    }
    return y;
  }
  
  
  real calc_log_prior(vector mu_vi,matrix Nu,matrix Ne,matrix Sp,matrix Sm,real c_v,real c_e,int num_r,int num_m,int num_ri,int num_mc,int[] idx_calc,int num_rc, vector v){
   matrix[num_rc,num_rc] Sigma_v;
   vector[num_rc] mu_v;
   real target_prior;
   mu_v = calc_mu(mu_vi,Nu,idx_calc,num_r);
   Sigma_v = calc_Sigma2(mu_vi,Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc);
   target_prior = multi_normal_lpdf(v | mu_v, Sigma_v);
   return target_prior;
 }
}

data{
    int<lower=0> num_r;
    int<lower=0> num_ri;
    int<lower=0> num_rd;
    int<lower=0> num_rc;
    int<lower=0> num_m;
    int<lower=0> num_mc;
    int<lower=0> num_b;
    int<lower=0> num_ri_wt;
    int<lower=0> num_smpl;
    matrix[num_m,num_r] S;
    matrix[num_r,num_ri] N;
    matrix[num_r,num_ri] Nu;
    matrix[num_r,num_rd] Ne;
    matrix[num_mc,num_r] Sp;
    matrix[num_mc,num_r] Sm;
    int idx_calc[num_rc];
    
    int<lower=0> num_g;
    int idx_g[num_g,2];

    int idx_b[num_b];
    
    vector[num_smpl] enz[num_rc];
    vector[num_smpl] sub1[num_rc];
    vector[num_smpl] pro1[num_rc];
    vector[num_smpl] sub2[num_rc];
    vector[num_smpl] pro2[num_rc];
    vector[num_smpl] met_eff[10];
    
    int is_irrev[num_rc];
    real<lower=0,upper=1> c_v;
    real<lower=0,upper=1> c_e;
    
    real v_max;
}
parameters{
    vector[num_rc] v[num_g];// flux
    real<lower=0> a[num_rc];
    real<lower=0> a_2;//TPS
    real<lower=0> b[num_b];
    real<lower=0> b_2;//FBA
    real<lower=0> e_atp[8];//ATP/ADP
    real<lower=0> e_nad[8];//NAD/NADH
    // real<lower=0,upper=1> e_nad[8];//NAD/NADH
    // real<lower=0,upper=1> e_nadp[4];//NADP/NADPH
    real<lower=0> e_allo[5];//HXK_T6P(i), PFK_F26bP(a),PFK_AMP(a),GPD_F16bP(i??),PYK_F16bP(a);
    vector<lower=-v_max,upper=v_max>[num_ri-1] mu_vi_wt_tmp;
    vector<lower=-v_max,upper=v_max>[num_ri] mu_vi_g[num_g-1];// mu_vi for other group
    real<lower=0> sigma_n;
    
    vector[num_rc] r_p[num_g];// protein turnover coefficient
    
    real<lower=0> sigma_p;// error term for prior in protein expression submodel
    
}

model{
  vector[num_smpl] v_x[num_rc];
  vector[num_ri] mu_vi_ref;
  vector[num_smpl] y[num_rc];
  vector[num_rc] mean_v;
  real b2[num_rc]  = rep_array(0.0,num_rc);
  vector[num_smpl] x[num_rc];
  
  b2[idx_b] = b;
      
  // calculate mean of all the prior mu_v
  mean_v = prep_mean_v(mu_vi_wt_tmp,mu_vi_g,
                           num_r,num_ri,num_rc,num_g,
                           Nu,idx_calc,idx_g);
  // calculate X
  x = calc_x(mean_v,sub1,sub2,pro1,pro2,met_eff,
                 a,a_2,b2,b_2,e_atp,e_nad,e_allo,
                 num_rc,num_smpl);                      
  
  for(r in 1:num_rc){
    for(g in 1:num_g){
      v_x[r][idx_g[g,1]:idx_g[g,2]] = v[g][r] * x[r][idx_g[g,1]:idx_g[g,2]];
    }
  }
  y = prep_y(r_p,v_x,idx_g,num_g,num_rc);
  
  for(r in 1:num_rc){
   target += normal_lpdf(enz[r] | v_x[r], sqrt(sigma_n) ); 
  }
  
  mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,num_ri);
  target += calc_log_prior(mu_vi_ref,Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[1]);
  for(g in 2:num_g){
    target += calc_log_prior(mu_vi_g[g-1],Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[g]);
  }
  
  target += normal_lpdf(a|0,1);
  target += normal_lpdf(a_2|0,1);
  target += normal_lpdf(b|0,1);
  target += normal_lpdf(b_2|0,1);
  target += normal_lpdf(e_atp|0,1);
  target += normal_lpdf(e_nad|0,1);
  target += normal_lpdf(e_allo|0,1);
  for(g in 1:num_g){
    target += normal_lpdf(r_p[g]|0,sigma_p);
  }
  target += cauchy_lpdf(sigma_p|0,0.5);
   
}
generated quantities{
  vector[num_smpl] log_lik[num_rc];
  real log_prior[num_g];
  //vector[num_smpl] enz_pred[num_rc];
  {
    vector[num_ri] mu_vi_ref;
    vector[num_rc] mean_v;
    real b2[num_rc]  = rep_array(0.0,num_rc);
    vector[num_smpl] x[num_rc];
    
    b2[idx_b] = b;
      
    // calculate mean of all the prior mu_v
    mean_v = prep_mean_v(mu_vi_wt_tmp,mu_vi_g,
                           num_r,num_ri,num_rc,num_g,
                           Nu,idx_calc,idx_g);
    // calculate X
    x = calc_x(mean_v,sub1,sub2,pro1,pro2,met_eff,
                 a,a_2,b2,b_2,e_atp,e_nad,e_allo,
                 num_rc,num_smpl);
    
    mu_vi_ref = prep_mu_vi_ref(mu_vi_wt_tmp,num_ri);
    
    log_prior[1] = calc_log_prior(mu_vi_ref,Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[1]);
    for(g in 2:num_g){
      log_prior[g] = calc_log_prior(mu_vi_g[g-1],Nu,Ne,Sp,Sm,c_v,c_e,num_r,num_m,num_ri,num_mc,idx_calc,num_rc,v[g]);
    }
    
    for(g in 1:num_g){
      for(r in 1:num_rc){
        for(n in idx_g[g,1]:idx_g[g,2]){
          log_lik[r][n] = normal_lpdf(enz[r][n] | v[g][r]*x[r][n], sqrt(sigma_n) );
          //enz_pred[r][n] = normal_rng(v[g][r]*x[r][n], sqrt(sigma_n));
        }
      }
    }
    
  }
}
