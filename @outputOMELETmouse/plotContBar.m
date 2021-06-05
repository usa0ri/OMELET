function plotContBar(obj,type,savedir)

    model_data = obj.model_data;
    mkdir(savedir);

    grp_names = model_data.grp_names;
    num_g = length(grp_names);
    cmb = obj.cont.cmb;

    opts = load_s_names(obj,type);

    for i=1:2
        is_errorbar = i==1;

        % inter-group
        for c=1:size(cmb,1)
            fname_str = [savedir '/' grp_names{cmb(c,1)} '_' grp_names{cmb(c,2)} ];
            plot_cont(opts.cont.intgrp(:,:,:,c),fname_str,model_data,is_errorbar,opts);
        end

    end
    
end

function plot_cont(cont_r_list,fname_str,model_data,is_errorbar,opts)

    fluxnames = model_data.X.rxn.rxn_names_include;
    num_rc = model_data.X.num.num_rc;
    num_r = length(opts.s_names);
    sz = size(cont_r_list);
    cont = nan(sz(1),num_r,sz(3));
    for i=1:num_r
       cont(:,i,:) = nansum(cont_r_list(:,opts.idx_s{i},:),2); 
    end
    cont_mu = nanmean(cont,3);
    cont_var = nanvar(cont,[],3);
    assert(all(sum(cont_mu,2)-1<1e-5));
    
    %%%%%%%%%%%%
    % reorder fluxnames
    % fluxnames = fluxnames(opts.idx_reorder);
    fluxnames = opts.fluxnames;
    cont_mu = cont_mu(opts.idx_reorder,:);
    cont_var = cont_var(opts.idx_reorder,:);
    %%%%%%%%%%%%
    
    % disp(cont_v_mat)
    fig  = figure('visible','off');
    % fig = figure();
    c_bar = categorical(fluxnames,fluxnames);
    bar_plot = bar(c_bar,cont_mu,'stacked','FaceColor','flat','EdgeColor','none');
    for i=1:num_r
        bar_plot(i).CData = opts.col(i,:);
    end
    hold on;
    if is_errorbar
        errorbar(cumsum(cont_mu')',sqrt(cont_var),'.k');
        ylim([0 1.5]);
    else
        ylim([0 1.0]);
    end
    
    legend(opts.s_names, 'Location', 'southeastoutside');

    fig.PaperUnits = 'inches';
    if is_errorbar
        fig.PaperSize = [num_rc 4];
        fig.PaperPosition = [0 0 num_rc 4];
%         saveas( gcf, [ fname_str '_error.png' ] );
        saveas( gcf, [ fname_str '_error.pdf' ] );
    else
        fig.PaperSize = [num_rc 2.5];
        fig.PaperPosition = [0 0 num_rc 2.5];
%         saveas( gcf, [ fname_str '.png' ] );
        saveas( gcf, [ fname_str '.pdf' ] );
    end
    close all;
end


function opts = load_s_names(obj,type)

    switch type
        case 'metab'
            cont = obj.cont;
           s_names =  {'Enzyme(protein)','Substrate','Product','Cofactor','Allosteric effectors','Unaccounted'};
           idx_s = {1,2,3,4:5,6:9,10};
           col = [lib.my_colors('matlab_blue',1,false);...
               lib.my_colors('matlab_orange',1,false);...
               lib.my_colors('matlab_yellow',1,false);...
               lib.my_colors('matlab_red',1,false);...
               lib.my_colors('matlab_purple',1,false);...
               [0.3 0.3 0.3];];
        case 'RNA'
            cont = obj.cont_rna;
            s_names = {'Transcripts','unaccounted'};
            idx_s = {1,2};
            col = [lib.my_colors('matlab_lblue',1,false);...
                [160 223 248]./255];
        case 'flux'
            cont = obj.cont_flux;
           s_names = {'Transcript','Unaccounted enzyme regulators',...
               'Substrate','Product','Cofactor','Allosteric effectors','Unaccounted flux regulators'}; 
           idx_s = {1,2,3,4,5:6,7:10,11};
           col = [lib.my_colors('matlab_lblue',1,false);...
               [160 223 248]./255;...
               lib.my_colors('matlab_orange',1,false);...
               lib.my_colors('matlab_yellow',1,false);...
               lib.my_colors('matlab_red',1,false);...
               lib.my_colors('matlab_purple',1,false);...
               [0.3 0.3 0.3];];

    end

    opts.idx_reorder = obj.fig_info.idx_reorder;
    opts.fluxnames = obj.fig_info.fluxnames;
    opts.cont = cont;
    opts.s_names = s_names;
    opts.idx_s = idx_s;
    opts.col = col;

end