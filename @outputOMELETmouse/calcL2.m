function calcL2(obj,savedir)
    
    mkdir(savedir);

    [data,v] = load_data(obj);
    
    calc_L2(data,obj,savedir);
    

end

function [data,v] = load_data(obj)
    
    model_data = obj.model_data;
    % reaction names
    rxn_names = model_data.X.rxn.rxn_names_include;
    num_rc = length(rxn_names);
    % reaction species
    idx_s = {1,2,3,4,5:6,7:10,11};
    s_names  = {'Enzyme(transcript)','Unaccounted enzyme regulator',...
        'Substrate','Product','Cofactor','Allosteric effectors',...
        'Unaccounted enzyme regualtor'}; 
    % s_names =  {'Et','Eu','S','P','C','A','U'};
    num_r = length(idx_s);
    % iterations
    iter = size(obj.par.a,1);
    % combinations of groups
    cmb = obj.cont.cmb;
    num_c = size(cmb,1);
    % group names
    grp_names = model_data.grp_names;
    num_g = length(grp_names);

    grp_names = model_data.grp_names;
    cmb_names = cell(1,num_c);
    for c=1:num_c
        grp1 = grp_names{cmb(c,1)};
        grp2 = grp_names{cmb(c,2)};
        cmb_names{c} = [grp1 '_' grp2];
    end

    % preparation for contribution data
    % cont_total = nan(num_rc,num_r,iter);
    cont_intgrp = nan(num_rc,num_r,iter,num_c);
    for c=1:num_c
        for i=1:num_r
            cont_intgrp(:,i,:,c) = nansum(obj.cont_flux.intgrp(:,idx_s{i},:,c),2);
        end
    end
   
    data.cont = cont_intgrp;
    data.dim1 = num_rc;
    data.dim1_names = rxn_names;
    data.dim2 = num_r;
    data.dim2_names = s_names;
    data.dim3 = iter;
    data.dim4 = num_c;
    data.dim4_names = cmb_names;
    data.grp_names = grp_names;
    data.num_grp = num_g;

% flux fold changes
    v = nan(iter,num_c,num_rc);
    for c=1:num_c
        v(:,c,:) = obj.par.v(:,cmb(c,2),:) ./...
            obj.par.v(:,cmb(c,1),:);
    end

end

function calc_L2(data,obj,savedir)
    
    x = data.cont;
    idx_cmb = [1 6;2 5];
    num_cmb = size(idx_cmb,1);
%     L1 = nan(data.dim1,data.dim3,num_cmb);
    L2 = nan(data.dim1,data.dim3,num_cmb);
    for c=1:num_cmb
%         L1(:,:,c) = nansum(abs(x(:,:,:,idx_cmb(c,1)) - x(:,:,:,idx_cmb(c,2))),2);
        L2(:,:,c) = sqrt(nansum((x(:,:,:,idx_cmb(c,1)) - x(:,:,:,idx_cmb(c,2))).^2,2));
    end
    
    plot_Lnorm(L2,data,obj,savedir);
    
end

function plot_Lnorm(L,data,obj,savedir)
    
    idx_cmb = [1 6;2 5];
    num_cmb = size(idx_cmb,1);

    % estimate distribution of L1-like distance
    % and compute 95% confidence(credit?) interval
    lowers = nan(data.dim1,num_cmb);
    uppers = nan(data.dim1,num_cmb);
    for c=1:num_cmb
        for r=1:data.dim1
            [f, xi] = ksdensity(L(r,:,c));
            f_ = f/sum(f);
            lower_now = max(xi(cumsum(f_)<0.01));
            upper_now = min(xi(cumsum(f_)>0.99));
            lowers(r,c) = lower_now;
            uppers(r,c) = upper_now;
        end
    end
    

    % xlabel is sorted
    idx_reorder = obj.fig_info.idx_reorder;
    L_reorder = L(idx_reorder,:,:);
    lowers_reorder = lowers(idx_reorder,:);
    uppers_reorder = uppers(idx_reorder,:);
    % violinplot
    for c=1:num_cmb
        fig = figure('visible','off');
        cmb_name = [data.dim4_names{idx_cmb(c,1)} '-' data.dim4_names{idx_cmb(c,2)}];
        lib.my_violin(L_reorder(:,:,c)',...
            'xlabel',obj.fig_info.fluxnames,...
            'facecolor',repmat([0.7 0.7 0.7],data.dim1,1),...
            'edgecolor',repmat([0.3 0.3 0.3],data.dim1,1),...
            'facealpha',1);
        ylim([0 sqrt(2)]);
        
        ylabel(['L2-like distance']);
        ax = gca;
        ax.XLim = ax.XLim + 0.025;
        ax.XTickLabelRotation = 30;
%         scatter(ax.XTick,lowers_reorder(:,c)',5,'+','k');
%         scatter(ax.XTick,uppers(:,c)','filled','k');
        
        title(cmb_name,'Interpreter','none');
        fig.PaperUnits = 'inches';
        fig.PaperSize = [5 1.5];
        fig.PaperPosition = [0 0 5 1.5];
        set(findobj(gca,'type','axes'),'FontSize',8,'FontName','SansSerif');
        print(gcf,'-painters',...
            [savedir '/violin_L2_' cmb_name '_reorder.pdf' ],...
            '-dpdf');
        close all;
    end
 

end