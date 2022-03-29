function calcContRNA2(obj)

    model_data = obj.model_data;
    par = obj.par;

    [dat_list,cont_rna] = calc_cont_mcmc(model_data,par);

    obj.cont_rna = cont_rna;
    obj.dat_cont_rna = dat_list;

end

function [dat_list, cont] = calc_cont_mcmc(model_data,par)

    %%% load data
    % dat: reactions * regulators * smpl * iter
    % dat_mu.total: reactions * regulators * iter

    % combination of groups for calculation of inter-grp contribution 
    num_g = model_data.out.num_g;
    cmb = nchoosek(1:num_g,2);

    [dat,dat_mu,dat_var] = load_data(model_data,par,cmb);

    [dat,dat_mu,dat_var] = calc_unknown(model_data,par,dat,dat_mu,dat_var,cmb);


    % pv_pr
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
    %     dat_var_intgrp_now = reshape(dat_var.intgrp(:,:,cmb(c,1),:) + dat_var.intgrp(:,:,cmb(c,2),:),...
    %         sz(1),sz(2),sz(3));
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
    num_rc = model_data.X.num.num_include;
    num_smpl = model_data.num_smpl;

    % sample index analyzed
    idx_tmp = sort(model_data.idx_g_org(:));
    idx_g = model_data.idx_g;
    num_g = model_data.num_g;

    % regulators
    % rna,unknown
    num_r=2;

    % rxn * reactants * smpl * iter
    dat = nan(num_rc,num_r,num_smpl,iter);

    %%% load mRNA
    data_rna = model_data.out.rna(:,idx_tmp);
    dat(:,1,:,:) = repmat(data_rna,1,1,iter);

    % calculate mu(total,grp) and variance(total,grp)
    % total : rxn * reactants * iter
    dat_mu_total = nan(num_rc,num_r,iter);
    dat_var_total = nan(num_rc,num_r,iter);
    dat_mu_grp = nan(num_rc,num_r,num_g,iter);
    dat_var_grp = nan(num_rc,num_r,num_g,iter);
    dat_var_intgrp = nan(num_rc,num_r,num_g,iter);
    for j=1:iter
        dat_now = reshape(dat(:,1,:,j),num_rc,num_smpl);
        [ mu_total_tmp, var_total_tmp ] = my_mu_var(dat_now);
        dat_mu_total(:,1,j) = mu_total_tmp;
        dat_var_total(:,1,j) = var_total_tmp;

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

        dat_mu_grp(:,1,:,j) = mu_g_tmp;
        dat_var_grp(:,1,:,j) = var_g_tmp; 
        dat_var_intgrp(:,1,:,j) = var_intg_tmp; 
    end

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
       dat_var_intgrp2(:,1,c,:) = (dat_mu.grp(:,1,cmb(c,1),:) - dat_mu.grp(:,1,cmb(c,2),:)).^2 ...
            .*((num_g1+num_g2)/num_smpl); 
    end

    dat_var.intgrp2 = dat_var_intgrp2;

end

function [dat,dat_mu,dat_var] = calc_unknown(model_data,par,dat,dat_mu,dat_var,cmb)

    iter = size(par.a,1);
    num_smpl = model_data.num_smpl;
    num_smpl_g = model_data.num_smpl_g;
    num_rc = model_data.X.num.num_include;
    num_g = model_data.num_g;
    idx_g = model_data.idx_g;

    r_p = 1+par.r_p_all;
    u = nan(iter,num_rc,num_smpl);
    for r=1:num_rc
       for g=1:num_g
    %       u(:,r,idx_g(g,1):idx_g(g,2)) =...
    %           par.enz_out(:,r,idx_g(g,1):idx_g(g,2)) -...
    %           nanmean(par.y(:,r,idx_g(g,1):idx_g(g,2)),[2 3]);      
            u(:,r,idx_g(g,1):idx_g(g,2)) =...
              par.enz_out(:,r,idx_g(g,1):idx_g(g,2)) -...
              reshape(model_data.out.rna(r,idx_g(g,1):idx_g(g,2)) .* r_p(:,g,r),...
              iter,1,num_smpl_g(g)); 
       end  
    end

    for j=1:iter
        dat_now = reshape(u(j,:,:),num_rc,num_smpl);
        dat(:,2,:,j) = dat_now;
        [ mu_total_tmp, var_total_tmp ] = my_mu_var(dat_now);
        dat_mu.total(:,2,j) = mu_total_tmp;
        dat_var.total(:,2,j) = var_total_tmp;

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

        dat_mu.grp(:,2,:,j) = mu_g_tmp;
        dat_var.grp(:,2,:,j) = var_g_tmp; 
        dat_var.intgrp(:,2,:,j) = var_intg_tmp;

    end

    % test
    for c=1:size(cmb,1)
        num_g1 = idx_g(cmb(c,1),2)-idx_g(cmb(c,1),1)+1;
        num_g2 = idx_g(cmb(c,2),2)-idx_g(cmb(c,2),1)+1;
       dat_var.intgrp2(:,2,c,:) = (dat_mu.grp(:,2,cmb(c,1),:) - dat_mu.grp(:,2,cmb(c,2),:)).^2 ...
            .*((num_g1+num_g2)/num_smpl); 
    end

end

function [mu_x,var_x] = my_mu_var(x)

    mu_x = nanmean(x,2);
    var_x = nansum((x-mu_x).^2,2)./size(x,2);

end

function pv_pr = prep_pvpr(model_data,par,dat_mu,cmb)

    iter = size(par.a,1);
    num_r = size(dat_mu.total,2);
    num_rc = size(dat_mu.total,1);
    num_g = model_data.out.num_g;

    num_c = size(cmb,1);
    for c=1:num_c
        dat_mu.cmb(:,:,c,:) = (dat_mu.grp(:,:,cmb(c,1),:) + dat_mu.grp(:,:,cmb(c,2),:))./2;
    end

    % sensitivity (total,grp,intgrp)
    pv_pr.total = nan(num_rc,num_r,iter);
    pv_pr.grp = nan(num_rc,num_r,iter,num_g);
    pv_pr.intgrp = nan(num_rc,num_r,iter,num_c);

    r_p = 1+par.r_p_all;
    for i=1:iter
        r_p_now = reshape(r_p(i,:,:),num_g,num_rc);

        pv_pr.total(:,1,i) = mean(r_p_now,1)';

        for g=1:num_g
            pv_pr.grp(:,1,i,g) = r_p_now(g,:)';
        end
        for c=1:num_c
            pv_pr.intgrp(:,1,i,c) = (r_p_now(cmb(c,1),:)+r_p_now(cmb(c,2),:))'./2;
        end
        % dE/dU
        pv_pr.total(:,2,i) = ones(num_rc,1,1);
        pv_pr.grp(:,2,i,:) = ones(num_rc,1,1,num_g);
        pv_pr.intgrp(:,2,i,:) = ones(num_rc,1,1,num_c);

    end

end

function cont_r = calc_cont(dat_var,pv_pr)

    num_rc = size(dat_var,1);
    num_r = size(dat_var,2);
    iter = size(dat_var,3);

    % calculate contribution of reactants
    % variance of v as sum of sensitivity*var of each reactants
    var_v_mat_sum = nan(iter,num_rc);
    % contribution of each reactants
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
        end
    end

end