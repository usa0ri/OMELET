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

    [dat,dat_mu,dat_var,met_eff_out] = load_data(model_data,par,cmb);

    [dat,dat_mu,dat_var] = calc_unknown(model_data,par,dat,dat_mu,dat_var,cmb);

    % pv_pr (sensitivity)
   pv_pr = prep_pvpr(model_data,par,dat_mu,cmb,met_eff_out);

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

function [dat,dat_mu,dat_var,met_eff_out] = load_data(model_data,par,cmb)

    iter = size(par.a,1);
    num_rc = model_data.X.num.num_include;
    num_smpl = model_data.num_smpl;
    r = model_data.X.rxn.rxn_names_include;
    m = model_data.X.met.met_names_eff;

    % sample index analyzed
    idx_tmp = sort(model_data.idx_g_org(:));
    idx_g = model_data.idx_g;
    num_g = model_data.num_g;
    num_smpl_g = idx_g(:,2)-idx_g(:,1)+1;

    % regulators:
    % enz,sub,pro
    regulators_name1 = {'enzyme','substrate','product'};
    if ~isempty(model_data.out.num_met_est)
        met_est_list = model_data.out.met_est_list;
        num_met_est_list = size(met_est_list,1);
        met_names_est = model_data.out.met_names_est;
        [~,idx_rxn_met_est] = ismember(met_est_list(:,2),r);
        [~,idx_met_met_est] = ismember(met_est_list(:,1),met_names_est);
        is_subpro_met_est = strcmp(met_est_list(:,3),'substrate');
    end

    % reactions * regulators * smpl * iter
    dat1 = nan(num_rc,length(regulators_name1),num_smpl,iter);

    %%% load enzyme(1), substrate(2), product(3)
    for i=1:iter
        dat1(:,1,:,i) = par.enz_out(i,:,:);
        dat1(:,2,:,i) = exp(model_data.out.sub(:,idx_tmp));
        dat1(:,3,:,i) = exp(model_data.out.pro(:,idx_tmp));
        for j=1:num_met_est_list
            for g=1:num_g-1
               dat1(idx_rxn_met_est(j),2,idx_g(g,1):idx_g(g,2),i) =...
                   repmat(par.met_est(i,idx_met_met_est(j),g),num_smpl_g(g),1); 
            end
            g=num_g;
            met_est_g = reshape(par.met_est(i,idx_met_met_est(j),:),num_g-1,1);
            met_est_now = (num_smpl-sum(met_est_g.*num_smpl_g(1:(g-1)),'all'))/num_smpl_g(g);
            assert(met_est_now>0);
            dat1(idx_rxn_met_est(j),2,idx_g(g,1):idx_g(g,2),i) =...
                   repmat(met_est_now,num_smpl_g(g),1); 
        end
    end
    assert(all(~isnan(dat1(:))));

    met_eff_list = model_data.D.met_eff_list_include;
    
    %%% load cofactors
    data_met_eff = exp(model_data.out.met_eff(:,idx_tmp));
    % stoichiometry matrix of cofactor
    met_eff_list_cofactor = met_eff_list(contains(met_eff_list(:,3),'cofactor'),:);
    met_eff_list_cofactor_sub = met_eff_list_cofactor(contains(met_eff_list_cofactor(:,3),'+'),:);
    met_eff_list_cofactor_pro = met_eff_list_cofactor(contains(met_eff_list_cofactor(:,3),'-'),:);
    
    S_eff = sparse(model_data.X.num.num_met_eff,model_data.X.num.num_include);
    S_cofactor_sub = S_eff;
    S_cofactor_pro = S_eff;
    for i=1:size(met_eff_list_cofactor_sub,1)
        S_cofactor_sub(ismember(m,met_eff_list_cofactor_sub(i,1)),...
            ismember(r,met_eff_list_cofactor_sub(i,2))) = 1;
    end
    for i=1:size(met_eff_list_cofactor_pro,1)
        S_cofactor_pro(ismember(m,met_eff_list_cofactor_pro(i,1)),...
            ismember(r,met_eff_list_cofactor_pro(i,2))) = 1;
    end
    
    met_eff_out = cell(2,3);
    met_eff_out(1,:) = {met_eff_list_cofactor_sub,S_cofactor_sub,'cofactor(sub)'};
    met_eff_out(2,:) = {met_eff_list_cofactor_pro,S_cofactor_pro,'cofactor(pro)'};
    
    %%% load allosteric effectors
    % stoichiometry matrix of allosteric effectors
    met_eff_list_allo = met_eff_list(contains(met_eff_list(:,3),'allosteric'),:);
    met_eff_list_allos = cell(2,1);
    met_eff_list_allos{1} = met_eff_list_allo(contains(met_eff_list_allo(:,3),'+'),:);
    met_eff_list_allos{2} = met_eff_list_allo(contains(met_eff_list_allo(:,3),'-'),:);
    
    % allosteric effectors: activators and inhibitors
    str_list = {'allosteric(a)','allosteric(i)'};
    for a=1:2
        str_now = str_list{a};
        met_eff_list_allo_now = met_eff_list_allos{a};
        cnt_rxn_a = cellfun(@(x) sum(ismember(met_eff_list_allo_now(:,2),x)),...
            r,'UniformOutput',true);
        num_allo_a = max(cnt_rxn_a);
        S_allo_a_list = cell(num_allo_a,1);
        for i=1:num_allo_a
           S_allo_a_list{i} = S_eff; 
        end
        met_allo_a_list = cell(num_allo_a,1);
        S_allo_a = S_eff;
        idx_cnt = [];
        for i=1:size(met_eff_list_allo_now,1)
            idx_rxn = ismember(r,met_eff_list_allo_now(i,2));
            cnt_rxn_a(idx_rxn) = cnt_rxn_a(idx_rxn)-1;
            if cnt_rxn_a(idx_rxn)>0
                cnt_now = cnt_rxn_a(idx_rxn)+1;
                S_allo_a_tmp = S_eff;
                S_allo_a_tmp(ismember(m,met_eff_list_allo_now(i,1)),...
                    ismember(r,met_eff_list_allo_now(i,2))) = 1;
                S_allo_a_list{cnt_now} = S_allo_a_list{cnt_now} + S_allo_a_tmp;
                met_allo_a_list{cnt_now} = [met_allo_a_list{cnt_now};met_eff_list_allo_now(i,:)];
                idx_cnt = [idx_cnt; i];
            else
                S_allo_a(ismember(m,met_eff_list_allo_now(i,1)),...
                    ismember(r,met_eff_list_allo_now(i,2))) = 1;
            end
        end
        S_allo_a_list{1} = S_allo_a;
        if isempty(idx_cnt)
            met_allo_a_list{1} = met_eff_list_allo_now;
        else
            met_allo_a_list{1} = met_eff_list_allo_now(all([1:size(met_eff_list_allo_now,1)]~=idx_cnt,1),:);
        end

        for i=1:num_allo_a
            met_eff_out = [met_eff_out;...
                {met_allo_a_list{i} S_allo_a_list{i},[str_now num2str(i)]}];
        end  
    end
    
    regulators_name2 = met_eff_out(:,3)';
    
    % reactions * regulators * smpl * iter
    dat2 = nan(num_rc,length(regulators_name2),num_smpl,iter);

    for i=1:iter
        for j=1:size(met_eff_out,1)
            dat2(:,j,:,i) = met_eff_out{j,2}' * data_met_eff;
        end
    end
    assert(all(~isnan(dat2(:))));

    
    dat = cat(2,dat1,dat2);

    regulators_name = [regulators_name1 regulators_name2];
    num_r = length(regulators_name);
    assert(num_r==size(dat,2));
    % calculate mu(total,grp) and variance(total,grp)
    % total : rxn * reactants * iter
    dat_mu_total = nan(num_rc,num_r,iter);
    dat_var_total = nan(num_rc,num_r,iter);
    dat_mu_grp = nan(num_rc,num_r,num_g,iter);
    dat_var_grp = nan(num_rc,num_r,num_g,iter);% within-group variance
    % between group variance: mu_grp - mu_total
    dat_var_intgrp = nan(num_rc,num_r,num_g,iter);
    for i=1:num_r
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

    assert(all(~isnan(dat_mu_total(:))));
    assert(all(~isnan(dat_var_grp(:))));
    assert(all(~isnan(dat_var_intgrp(:))));

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
       dat_var_intgrp2(:,:,c,:) = (dat_mu.grp(:,:,cmb(c,1),:) - dat_mu.grp(:,:,cmb(c,2),:)).^2 ...
            .*((num_g1+num_g2)/num_smpl); 
    end
    dat_var.intgrp2 = dat_var_intgrp2;

end

function [dat,dat_mu,dat_var] = calc_unknown(model_data,par,dat,dat_mu,dat_var,cmb)

    iter = size(par.a,1);
    mean_v = par.mean_v';

    dat_mu.v = mean_v;

    num_smpl = model_data.num_smpl;
    num_rc = model_data.X.num.num_include;
    num_g = model_data.num_g;
    idx_g = model_data.idx_g;

    v_k = par.enz_pred ./ par.x;
    u = nan(num_rc,1,num_smpl,iter);
    for g=1:num_g
        v_now = reshape(par.v(:,g,:),iter,num_rc);
        u_tmp = v_k(:,:,idx_g(g,1):idx_g(g,2)) - ...
           repmat(v_now,1,1,(idx_g(g,2)-idx_g(g,1)+1));
       for i=idx_g(g,1):idx_g(g,2)
           u(:,1,i,:) = u_tmp(:,:,i-idx_g(g,1)+1)';
       end
    end

    dat = cat(2,dat,u);

    dat_mu_total = nan(num_rc,1,iter);
    dat_var_total = nan(num_rc,1,iter);
    dat_mu_grp = nan(num_rc,1,num_g,iter);
    dat_var_grp = nan(num_rc,1,num_g,iter);
    dat_var_intgrp = nan(num_rc,1,num_g,iter);
    for i=1:iter
        [ mu_total_tmp, var_total_tmp ] = my_mu_var(reshape(u(:,1,:,i),num_rc,num_smpl));

        dat_mu_total(:,1,i) = mu_total_tmp;
        dat_var_total(:,1,i) = var_total_tmp;

        mu_g_tmp = nan(num_rc,num_g);
        var_g_tmp = nan(num_rc,num_g);
        var_intg_tmp = nan(num_rc,num_g);
        for g=1:num_g
            dat_g_ = u(:,1,idx_g(g,1):idx_g(g,2),i);
            [ mu_g_tmp_, var_g_tmp_ ] = my_mu_var(reshape(dat_g_,num_rc,size(dat_g_,3)));
            mu_g_tmp(:,g) = mu_g_tmp_;
            var_g_tmp(:,g) = var_g_tmp_;
            var_intg_tmp(:,g) = (mu_total_tmp - mu_g_tmp_).^2 .*((idx_g(g,2)-idx_g(g,1)+1)/num_smpl);
        end

        dat_mu_grp(:,1,:,i) = mu_g_tmp;
        dat_var_grp(:,1,:,i) = var_g_tmp; 
        dat_var_intgrp(:,1,:,i) = var_intg_tmp; 

    end
    
    dat_mu.total = cat(2,dat_mu.total,dat_mu_total);
    dat_var.total = cat(2,dat_var.total,dat_var_total);
    dat_mu.grp = cat(2,dat_mu.grp,dat_mu_grp);
    dat_var.grp = cat(2,dat_var.grp,dat_var_grp);
    dat_var.intgrp = cat(2,dat_var.intgrp,dat_var_intgrp);

    assert(~any(isnan(dat(:))));

    dat_var_intgrp2 = nan(num_rc,1,num_g,iter);
    for c=1:size(cmb,1)
        num_g1 = idx_g(cmb(c,1),2)-idx_g(cmb(c,1),1)+1;
        num_g2 = idx_g(cmb(c,2),2)-idx_g(cmb(c,2),1)+1;
       dat_var_intgrp2(:,1,c,:) = (dat_mu_grp(:,1,cmb(c,1),:) - dat_mu_grp(:,1,cmb(c,2),:)).^2 ...
            .*((num_g1+num_g2)/num_smpl); 
    end
    dat_var.intgrp2 = cat(2,dat_var.intgrp2,dat_var_intgrp2);

end

function [mu_x,var_x] = my_mu_var(x)

    mu_x = nanmean(x,2);
    var_x = nansum((x-mu_x).^2,2)./size(x,2);

end

function pv_pr = prep_pvpr(model_data,par,dat_mu,cmb,met_eff_out)

    iter = size(par.a,1);
    r = model_data.X.rxn.rxn_names_include;
    num_r = size(dat_mu.total,2);
    num_rc = size(dat_mu.total,1);
    num_g = model_data.num_g;
    
    el = zeros(num_rc,num_r-2,iter);
    
    % elasticity for substrates
    el(:,1,:) = par.a';
    
    % elasticity for products
    el(:,2,:) = par.b';
    
    % elasticity for cofactors and allosteric effectors
    met_eff_list = model_data.out.met_eff_list_include;
    idx_total = false(size(met_eff_list,1),1);
    for i=1:size(met_eff_out,1)
       met_eff_now = met_eff_out{i,1};
       for j=1:size(met_eff_now,1)
           idx_now = ismember(met_eff_list(:,2),met_eff_now{j,2}) &...
               ismember(met_eff_list(:,1),met_eff_now{j,1});
           sign_now = met_eff_list{j,3}(1);
           if strcmp(sign_now,'+')
               el(ismember(r,met_eff_now{j,2}),i+2,:) = par.e_effs(:,idx_now);
           elseif strcmp(sign_now,'-')
               el(ismember(r,met_eff_now{j,2}),i+2,:) = -par.e_effs(:,idx_now);
           end
           idx_total(idx_now) = true;
       end
    end
    assert(all(idx_total));
    
    % mu of group combinations (cmb)
    num_c = size(cmb,1);
    for c=1:num_c
        dat_mu.cmb(:,:,c,:) = (dat_mu.grp(:,:,cmb(c,1),:) + dat_mu.grp(:,:,cmb(c,2),:))./2;
    end
    
    
    % sensitivity (total,grp,intgrp)
    pv_pr.total = nan(num_rc,num_r,iter);
    pv_pr.grp = nan(num_rc,num_r,iter,num_g);
    pv_pr.intgrp = nan(num_rc,num_r,iter,num_c);
    idx_met = 2:(num_r-1);
    for i=1:iter
        % dv/dE
        pv_pr.total(:,1,i) = dat_mu.v(:,i) .* (1 + nansum(el(:,:,i).*log(dat_mu.total(:,idx_met,i)),2));
        for g=1:num_g
            pv_pr.grp(:,1,i,g) = dat_mu.v(:,i) .* (1 + nansum(el(:,:,i).*log(dat_mu.grp(:,idx_met,g,i)),2));
        end
        for c=1:num_c
            pv_pr.intgrp(:,1,i,c) = dat_mu.v(:,i) .* (1 + nansum(el(:,:,i).*log(dat_mu.cmb(:,idx_met,c,i)),2));
        end
        
        % dv/dM
        for j=1:(num_r-2)
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
        pv_pr.total(:,num_r,i) = ones(num_rc,1,1);
        pv_pr.grp(:,num_r,i,:) = ones(num_rc,1,1,num_g);
        pv_pr.intgrp(:,num_r,i,:) = ones(num_rc,1,1,num_c);
        
    end
    pv_pr.total(isnan(pv_pr.total)) = 0;
    pv_pr.total(isinf(pv_pr.total)) = 0;
    
    assert(all(~isnan(pv_pr.total(:))));

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
