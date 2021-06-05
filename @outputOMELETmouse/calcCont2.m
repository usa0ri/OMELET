function calcCont2(obj)

    model_data = obj.model_data;
    par = obj.par;

    [dat,cont] = calc_cont_mcmc(model_data,par);

    obj.cont = cont;
    obj.dat_cont = dat;

end

function [dat_list,cont] = calc_cont_mcmc(model_data,par)

    %%% load data
    % dat: reactions * regulators * smpl * iter
    % dat_mu.total: reactions * regulators * iter

    % combination of groups for calculation of inter-grp contribution 
    num_g = model_data.num_g;
    cmb = nchoosek(1:num_g,2);

    [dat,dat_mu,dat_var] = load_data(model_data,par,cmb);

    [dat,dat_mu,dat_var] = calc_unknown(model_data,par,dat,dat_mu,dat_var,cmb);

    % pv_pr (sensitivity)
    pv_pr = prep_pvpr(model_data,par,dat_mu,cmb);

    % total
    cont.total = calc_cont(dat_var.total,pv_pr.total);

    % each grp
    sz = size(dat_var.total);
    for g = 1:num_g
        dat_var_grp_now = reshape(dat_var.grp(:,:,g,:),sz(1),sz(2),sz(3));
        cont.grp(:,:,:,g) = calc_cont(dat_var_grp_now,pv_pr.grp(:,:,:,g));
    end

    % inter-grp
    cont.cmb = cmb;
    for c=1:size(cmb,1)
        dat_var_intgrp_now = reshape(dat_var.intgrp2(:,:,c,:),...
            sz(1),sz(2),sz(3));
        cont.intgrp(:,:,:,c) = calc_cont(dat_var_intgrp_now,pv_pr.intgrp(:,:,:,c));
    end

    dat_list.dat = dat;
    dat_list.dat_mu = dat_mu;
    dat_list.dat_var = dat_var;
    dat_list.pv_pr = pv_pr;

end

function [dat,dat_mu,dat_var] = load_data(model_data,par,cmb)

    iter = size(par.a,1);
    num_rc = model_data.X.num.num_rc;
    num_smpl = model_data.num_smpl;
    r = model_data.X.rxn.rxn_names_include;
    m = model_data.X.met.met_names_eff;

    % sample index analyzed
    idx_tmp = sort(model_data.idx_g_org(:));
    idx_g = model_data.idx_g;
    num_g = model_data.num_g;

    % regulators:
    % enz,sub,pro,co_sub,co_pro,allo_a,allo_i1,allo_i2,allo_i3,unknown
    num_r=10;

    % reactions * regulators * smpl * iter
    dat = nan(num_rc,num_r,num_smpl,iter);

    %%% load enzyme(1), substrate(2), product(3)
    for i=1:iter
        dat(:,1,:,i) = par.enz_out(i,:,:);
        dat(:,2,:,i) = exp(model_data.out.sub1(:,idx_tmp));
        dat(:,3,:,i) = exp(model_data.out.pro1(:,idx_tmp));
    %     [~,idx_OAA_sub] = ismember({'Pck1','Cs'},r);
        [~,idx_OAA_sub] = ismember('Pck1',r);
        [~,idx_OAA_pro] = ismember({'Mdh2','Cs'},r); % Cs has 2 substrates
        [~,idx_Pyr_sub] = ismember('Pcx',r);
        [~,idx_Pyr_pro] = ismember({'Ldha','Gpt'},r);
        for g=1:num_g
            % OAA_Pck1 (sub)
            dat(idx_OAA_sub,2,idx_g(g,1):idx_g(g,2),i) = par.c_OAA_out(i,g);
            % OAA_Mdh2 (pro), OAA_Cs (sub)
            dat(idx_OAA_pro,3,idx_g(g,1):idx_g(g,2),i) = par.c_OAA_out(i,g);
            % Pyr_Pcx (sub)
            dat(idx_Pyr_sub,2,idx_g(g,1):idx_g(g,2),i) = par.c_Pyr_out(i,g);
            % Pyr_Ldha (pro), Pyr_Gpt (pro)
            dat(idx_Pyr_pro,3,idx_g(g,1):idx_g(g,2),i) = par.c_Pyr_out(i,g);
        end
    end


    %%% load cofactor(4,5)
    data_met_eff = exp(model_data.out.met_eff(:,idx_tmp));
    % stoichiometry matrix of cofactor
    S_eff = sparse(model_data.X.num.num_met_eff,model_data.X.num.num_rc);
    S_cofactor_sub = S_eff;
    S_cofactor_pro = S_eff;
    S_cofactor_sub(ismember(m,'NAD+'),ismember(r,'Gpd1')) = 1;
    S_cofactor_pro(ismember(m,'NADH'),ismember(r,'Gpd1')) = 1;
    S_cofactor_sub(ismember(m,'NAD+'),ismember(r,'Ldha')) = 1;
    S_cofactor_pro(ismember(m,'NADH'),ismember(r,'Ldha')) = 1;
    S_cofactor_sub(ismember(m,'NAD+'),ismember(r,'Mdh2')) = 1;
    S_cofactor_pro(ismember(m,'NADH'),ismember(r,'Mdh2')) = 1;
    S_cofactor_sub(ismember(m,'NAD+'),ismember(r,'Glud1')) = 1;
    S_cofactor_pro(ismember(m,'NADH'),ismember(r,'Glud1')) = 1;
    S_cofactor_sub(ismember(m,'ADP'),ismember(r,'Pklr')) = 1;
    S_cofactor_sub(ismember(m,'ATP'),ismember(r,'Pcx')) = 1;
    S_cofactor_sub(ismember(m,'GTP'),ismember(r,'Pck1')) = 1;
    S_cofactor_sub(ismember(m,'FAD'),ismember(r,'Sdha')) = 1;


    %%% load allosteric effectors (6,7)
    % stoichiometry matrix of allosteric effectors
    S_allo_i = S_eff;
    S_allo_a = S_eff;
    S_allo_a(ismember(m,'Cit'),ismember(r,'Fbp1')) = 1;
    S_allo_i(ismember(m,'AMP'),ismember(r,'Fbp1')) = 1;
    S_allo_a(ismember(m,'F1,6P'),ismember(r,'Pklr')) = 1;
    S_allo_i(ismember(m,'ATP'),ismember(r,'Pklr')) = 1;
    S_allo_a(ismember(m,'AcCoA'),ismember(r,'Pcx')) = 1;
    S_allo_a(ismember(m,'ADP'),ismember(r,'Glud1')) = 1;
    S_allo_i(ismember(m,'ATP'),ismember(r,'Glud1')) = 1;

    S_allo_i2 = S_eff;
    S_allo_i2(ismember(m,'Ala'),ismember(r,'Pklr')) = 1;
    S_allo_i2(ismember(m,'GTP'),ismember(r,'Glud1')) = 1;

    S_allo_i3 = S_eff;
    S_allo_i3(ismember(m,'Phe'),ismember(r,'Pklr')) = 1;
    S_allo_i3(ismember(m,'Leu'),ismember(r,'Glud1')) = 1;


    for i=1:iter
        dat(:,4,:,i) = S_cofactor_sub' * data_met_eff;
        dat(:,5,:,i) = S_cofactor_pro' * data_met_eff;
        dat(:,6,:,i) = S_allo_a' * data_met_eff;
        dat(:,7,:,i) = S_allo_i' * data_met_eff;
        dat(:,8,:,i) = S_allo_i2' * data_met_eff;
        dat(:,9,:,i) = S_allo_i3' * data_met_eff;
    end

    % calculate mu(total,grp) and variance(total,grp)
    % total : rxn * reactants * iter
    dat_mu_total = nan(num_rc,num_r,iter);
    dat_var_total = nan(num_rc,num_r,iter);
    dat_mu_grp = nan(num_rc,num_r,num_g,iter);
    dat_var_grp = nan(num_rc,num_r,num_g,iter);% within-group variance
    % between group variance: mu_grp - mu_total
    dat_var_intgrp = nan(num_rc,num_r,num_g,iter);
    for i=1:num_r-1% except unknown term
        for j=1:iter
            dat_now = reshape(dat(:,i,:,j),num_rc,num_smpl);
            [ mu_total_tmp, var_total_tmp ] = my_mu_var(dat_now);
            dat_mu_total(:,i,j) = mu_total_tmp;
            dat_var_total(:,i,j) = var_total_tmp;

            mu_g_tmp = nan(num_rc,num_g);
            var_g_tmp = nan(num_rc,num_g);
            var_intg_tmp = nan(num_rc,num_g);
            for g=1:num_g
                dat_g_ = reshape(dat_now(:,idx_g(g,1):idx_g(g,2)),num_rc,(idx_g(g,2)-idx_g(g,1)+1));
                [ mu_g_tmp_, var_g_tmp_ ] = my_mu_var(dat_g_);
                mu_g_tmp(:,g) = mu_g_tmp_;
                var_g_tmp(:,g) = var_g_tmp_;
                var_intg_tmp(:,g) = (mu_total_tmp - mu_g_tmp_).^2 .*((idx_g(g,2)-idx_g(g,1)+1)/num_smpl);
            end

            dat_mu_grp(:,i,:,j) = mu_g_tmp;
            dat_var_grp(:,i,:,j) = var_g_tmp; 
            dat_var_intgrp(:,i,:,j) = var_intg_tmp; 
        end
    end

    % except unknown(10)
    tmp = dat_mu_total(:,1:(num_r-1),:,:);
    assert(all(~isnan(tmp(:))));
    tmp = dat_var_grp(:,1:(num_r-1),:,:);
    assert(all(~isnan(tmp(:))));
    tmp = dat_var_intgrp(:,1:(num_r-1),:,:);
    assert(all(~isnan(tmp(:))));

    dat_mu.total = dat_mu_total;
    dat_mu.grp = dat_mu_grp;
    dat_var.total = dat_var_total;
    dat_var.grp = dat_var_grp;
    dat_var.intgrp = dat_var_intgrp;

    % between-group variance is evaluated by mu_grp1 - mu_grp2
    dat_var_intgrp2 = nan(num_rc,num_r,size(cmb,1),iter);
    for c=1:size(cmb,1)
        num_g1 = idx_g(cmb(c,1),2)-idx_g(cmb(c,1),1)+1;
        num_g2 = idx_g(cmb(c,2),2)-idx_g(cmb(c,2),1)+1;
       dat_var_intgrp2(:,1:9,c,:) = (dat_mu.grp(:,1:9,cmb(c,1),:) - dat_mu.grp(:,1:9,cmb(c,2),:)).^2 ...
            .*((num_g1+num_g2)/num_smpl); 
    end
    dat_var.intgrp2 = dat_var_intgrp2;

end

function [dat,dat_mu,dat_var] = calc_unknown(model_data,par,dat,dat_mu,dat_var,cmb)

    iter = size(par.a,1);
    mean_v = par.mean_v';

    dat_mu.v = mean_v;

    num_smpl = model_data.num_smpl;
    num_rc = model_data.X.num.num_rc;
    num_g = model_data.num_g;
    idx_g = model_data.idx_g;

    v_k = par.enz_out ./ par.x;
    u = nan(num_rc,num_smpl,iter);
    for g=1:num_g
        v_now = reshape(par.v(:,g,:),iter,num_rc);
        u_tmp = v_k(:,:,idx_g(g,1):idx_g(g,2)) - ...
           repmat(v_now,1,1,(idx_g(g,2)-idx_g(g,1)+1));
       for i=idx_g(g,1):idx_g(g,2)
           u(:,i,:) = u_tmp(:,:,i-idx_g(g,1)+1)';
       end
    end

    dat(:,10,:,:) = u;

    for i=1:iter
        [ mu_total_tmp, var_total_tmp ] = my_mu_var(u(:,:,i));

        dat_mu.total(:,10,i) = mu_total_tmp;
        dat_var.total(:,10,i) = var_total_tmp;

        mu_g_tmp = nan(num_rc,num_g);
        var_g_tmp = nan(num_rc,num_g);
        var_intg_tmp = nan(num_rc,num_g);
        for g=1:num_g
            dat_g_ = u(:,idx_g(g,1):idx_g(g,2),i);
            [ mu_g_tmp_, var_g_tmp_ ] = my_mu_var(dat_g_);
            mu_g_tmp(:,g) = mu_g_tmp_;
            var_g_tmp(:,g) = var_g_tmp_;
            var_intg_tmp(:,g) = (mu_total_tmp - mu_g_tmp_).^2 .*((idx_g(g,2)-idx_g(g,1)+1)/num_smpl);
        end

        dat_mu.grp(:,10,:,i) = mu_g_tmp;
        dat_var.grp(:,10,:,i) = var_g_tmp; 
        dat_var.intgrp(:,10,:,i) = var_intg_tmp; 

    end

    assert(~any(isnan(dat(:))));

    for c=1:size(cmb,1)
        num_g1 = idx_g(cmb(c,1),2)-idx_g(cmb(c,1),1)+1;
        num_g2 = idx_g(cmb(c,2),2)-idx_g(cmb(c,2),1)+1;
       dat_var.intgrp2(:,10,c,:) = (dat_mu.grp(:,10,cmb(c,1),:) - dat_mu.grp(:,10,cmb(c,2),:)).^2 ...
            .*((num_g1+num_g2)/num_smpl); 
    end

end

function [mu_x,var_x] = my_mu_var(x)

    mu_x = nanmean(x,2);
    var_x = nansum((x-mu_x).^2,2)./size(x,2);

end

function pv_pr = prep_pvpr(model_data,par,dat_mu,cmb)

    iter = size(par.a,1);
    r = model_data.X.rxn.rxn_names_include;
    num_r = size(dat_mu.total,2);
    num_rc = size(dat_mu.total,1);
    num_g = model_data.num_g;
    
    el = zeros(num_rc,8,iter);
    
    % elasticity for substrates
    el(:,1,:) = par.a';
    
    % elasticity for products
    el(:,2,:) = par.b';
    
    % elasticity for cofactor (3-4)
    el(ismember(r,'Gpd1'),3,:) = par.e_cofactor(:,1);
    el(ismember(r,'Gpd1'),4,:) = -par.e_cofactor(:,2);
    el(ismember(r,'Ldha'),3,:) = par.e_cofactor(:,3);
    el(ismember(r,'Ldha'),4,:) = -par.e_cofactor(:,4);
    el(ismember(r,'Mdh2'),3,:) = par.e_cofactor(:,5);
    el(ismember(r,'Mdh2'),4,:) = -par.e_cofactor(:,6);
    el(ismember(r,'Glud1'),3,:) = par.e_cofactor(:,7);
    el(ismember(r,'Glud1'),4,:) = -par.e_cofactor(:,8);
    el(ismember(r,'Pklr'),3,:) = par.e_cofactor(:,9);
    % el(ismember(r,'Pklr'),4,:) = -par.e_cofactor(:,10); -> allosteric
    el(ismember(r,'Pcx'),3,:) = par.e_cofactor(:,11);
    el(ismember(r,'Pck1'),3,:) = par.e_cofactor(:,12);
    el(ismember(r,'Sdha'),3,:) = par.e_cofactor(:,13);
    
    % elasticity for allostery (5-8)
    el(ismember(r,'Fbp1'),5,:) = par.e_allo(:,2);
    el(ismember(r,'Fbp1'),6,:) = -par.e_allo(:,1);
    el(ismember(r,'Pklr'),5,:) = par.e_allo(:,3);
    el(ismember(r,'Pklr'),6,:) = -par.e_cofactor(:,10);
    el(ismember(r,'Pcx'),5,:) = par.e_allo(:,6);
    el(ismember(r,'Glud1'),5,:) = par.e_allo(:,8);
    el(ismember(r,'Glud1'),6,:) = -par.e_allo(:,7);
    
    el(ismember(r,'Pklr'),7,:) = par.e_allo(:,4);
    el(ismember(r,'Glud1'),7,:) = -par.e_allo(:,9);
    
    el(ismember(r,'Pklr'),8,:) = -par.e_allo(:,5);
    el(ismember(r,'Glud1'),8,:) = -par.e_allo(:,10);
    
    
    % mu of group combinations (cmb)
    num_c = size(cmb,1);
    for c=1:num_c
        dat_mu.cmb(:,:,c,:) = (dat_mu.grp(:,:,cmb(c,1),:) + dat_mu.grp(:,:,cmb(c,2),:))./2;
    end
    
    
    % sensitivity (total,grp,intgrp)
    pv_pr.total = nan(num_rc,num_r,iter);
    pv_pr.grp = nan(num_rc,num_r,iter,num_g);
    pv_pr.intgrp = nan(num_rc,num_r,iter,num_c);
    for i=1:iter
        % dv/dE
        pv_pr.total(:,1,i) = dat_mu.v(:,i) .* (1 + nansum(el(:,:,i).*log(dat_mu.total(:,2:9,i)),2));
        for g=1:num_g
            pv_pr.grp(:,1,i,g) = dat_mu.v(:,i) .* (1 + nansum(el(:,:,i).*log(dat_mu.grp(:,2:9,g,i)),2));
        end
        for c=1:num_c
            pv_pr.intgrp(:,1,i,c) = dat_mu.v(:,i) .* (1 + nansum(el(:,:,i).*log(dat_mu.cmb(:,2:9,c,i)),2));
        end
        
        % dv/dM
        for j=1:8
            pv_pr.total(:,j+1,i) = ( dat_mu.v(:,i) .* dat_mu.total(:,1,i) .* el(:,j,i) )...
                ./ dat_mu.total(:,j+1,i);
            for g=1:num_g
                pv_pr.grp(:,j+1,i,g) = ( dat_mu.v(:,i) .* dat_mu.grp(:,1,g,i) .* el(:,j,i) )...
                    ./ dat_mu.grp(:,j+1,g,i);
            end
            for c=1:num_c
                pv_pr.intgrp(:,j+1,i,c) = ( dat_mu.v(:,i) .* dat_mu.cmb(:,1,c,i) .* el(:,j,i) )...
                    ./ dat_mu.cmb(:,j+1,c,i);
            end
        end
        % dv/dU
        pv_pr.total(:,10,i) = ones(num_rc,1,1);
        pv_pr.grp(:,10,i,:) = ones(num_rc,1,1,num_g);
        pv_pr.intgrp(:,10,i,:) = ones(num_rc,1,1,num_c);
        
    end

end

function cont_r = calc_cont(dat_var,pv_pr)

    num_rc = size(dat_var,1);
    num_r = size(dat_var,2);
    iter = size(dat_var,3);

    % calculate contribution of regulators
    % variance of v as sum of sensitivity*var of each regulator
    var_v_mat_sum = nan(iter,num_rc);
    % contribution of each regulators
    cont_r = nan(size(pv_pr));
    for s = 1:iter
        var_r_now = dat_var(:,:,s);
        pv_pr_now = pv_pr(:,:,s);
        pv_pr_now(isnan(pv_pr_now)) = 0;
        for i=1:num_rc
            Sigma_now = diag(var_r_now(i,:));
            var_v_mat_sum(i,s) = pv_pr_now(i,:) * Sigma_now * pv_pr_now(i,:)';
            for j=1:num_r
                cont_r(i,j,s) = ( pv_pr_now(i,j) * Sigma_now(j,j) * pv_pr_now(i,j) )...
                    ./ var_v_mat_sum(i,s);
            end
            if i==12% Cs
                % combine contributions of substrate and products in Cs
                % contributions of substrate: AcCoA
                % contributions of product: OAA
                cont_r(i,2,s) = cont_r(i,2,s) + cont_r(i,3,s);% substrate
                cont_r(i,3,s) = 0;% product

            end
        end
    end


end
