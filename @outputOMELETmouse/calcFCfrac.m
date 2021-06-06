function calcFCfrac(obj,savedir)

    mkdir(savedir);

    % G6pc*6 = Pgm2*6 + Gpd1*3 + Eno1*3
    % rxn_frac1 = {'Pgm2','Gpd1','Eno1'};
    rxn_frac1 = {'Eno1','Gpd1','Pgm2'};

    % Pck1*3 = Pklr*3 + Eno1*3
    rxn_frac2 = {'Pklr','Eno1'};

    % Pcx*3 = Pklr*3 + Ldha*3 + Gpt*3
    rxn_frac3 = {'Pklr','Ldha','Gpt'};

    % G6pc*6 = Pgm2*6 + Gpd1*3 + Eno1*3
    %        =                   Pck1*3 - Pklr*3
    %        =                   Pcx*3 + Mdh2* - Cs* - Pklr*3
    %        =                   Glud1*3 + Ldha*3 + Gpt*3
    % rxn_frac4 = {'Pgm2','Gpd1','Ldha','Gpt','Glud1'};
    rxn_frac4 = {'Glud1','Gpt','Ldha','Gpd1','Pgm2'};

    rxn_frac_list = {rxn_frac1, rxn_frac2, rxn_frac3, rxn_frac4};
    rxn_target_list = {'HGPold','Pck1','Pcx','HGP'};
    for i=1:length(rxn_target_list)
        plot_frac_bar(rxn_frac_list{i},rxn_target_list{i},obj,savedir);
    end

end

function plot_frac_bar(rxn_frac,rxn_target,obj,savedir)

    grp_names = obj.model_data.grp_names;
    num_g = length(grp_names);
    fluxnames = obj.model_data.X.rxn.rxn_names_include;
    [~,idx_frac] = ismember(rxn_frac,fluxnames);
    [~,idx_target] = ismember(rxn_target,fluxnames);

    v = obj.par.v(:,:,idx_frac);

    if contains(rxn_target,'HGP')
        v(:,:,end) = 2.*v(:,:,end);
        vt = sum(v,3);
    else
    %     vt = obj.par.v(:,:,idx_target);
        vt = sum(v,3);
    end
    frac = nan(size(v));
    for i=1:size(v,3)
        frac(:,:,i) = v(:,:,i)./vt;
    end

    stat_names = {'Median','Mean','SD',...
        'Lower limit of 95% credible interval',...
        'Upper limit of 95% credible interval'};
    num_stat = length(stat_names);
    frac_median = median(frac,1);
    frac_mean = mean(frac,1);
    frac_std = std(frac,0,1);
    frac_95lower = prctile(frac,2.5,1);
    frac_95upper = prctile(frac,97.5,1);
    frac_stat = cat(1,frac_median,frac_mean,frac_std,frac_95lower,frac_95upper);

    colname_frac = cellfun(@(y) cellfun(@(x) [y ' (' x ')'], stat_names,'UniformOutput',false),...
        grp_names,'UniformOutput',false);
    rowname_frac = repmat([{''} rxn_frac]',num_g,1);

    % reshape data_now
    data_tmp = {};
    for g=1:num_g
       data_tmp = [data_tmp; colname_frac{g};...
           num2cell(reshape(frac_stat(:,g,:),num_stat,length(idx_frac))')];
    end
    out = [[{''}; rowname_frac], [ stat_names; data_tmp]];
    writecell(out,[savedir '/data_' rxn_target '.csv']);

    % plot
    frac_mu = reshape(mean(frac,1),num_g,length(idx_frac));
    frac_std = reshape(frac_std,num_g,length(idx_frac));

    fig = figure('visible','off');
    c = categorical(grp_names,grp_names);
    bar(c,frac_mu,...
        'stacked','FaceColor','flat');
    b_tmp = findobj(gcf,'type','bar');
    if length(idx_frac)==2
        col_list = {[0.8 0.8 0.8],[0.3 0.3 0.3]};
    elseif length(idx_frac)==3
        col_list = {[0.8 0.8 0.8],[0.5 0.5 0.5],[0.3 0.3 0.3]};
    elseif length(idx_frac)==5
        col_list = {[1 1 1],[0.8 0.8 0.8],[0.5 0.5 0.5],[0.3 0.3 0.3],[0.1 0.1 0.1]};
    end
    for i=1:length(b_tmp)
        b_tmp(i).FaceColor = col_list{i};
    end
    hold on;
    % flip legends
    bb = fig.Children.Children';
    l = legend(bb,rxn_frac,'Location','southeastoutside');
    legholder = l.String;
    for i=1:length(legholder)
       l.String(i) = legholder(end+1-i);
    end

    errorbar(cumsum(frac_mu')',frac_std,'.k');

    ylabel(['Fractional contribution to ' rxn_target ]);

    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [5 3];
    fig.PaperPosition = [0 0 5 3];
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');

    % saveas( gcf, [ savedir '/v_frac' rxn_target '.png' ] );
    saveas( gcf, [ savedir '/v_frac' rxn_target '.pdf' ] );
    close all;

end
