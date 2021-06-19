function cmpTracer(obj,savedir)

    mkdir(savedir);

    out = load_data;
    plot_tracer(obj,out,savedir);
    % plot_tracer_abs(obj,out,savedir);

    scatter_tracer(obj,out,savedir);

    out_ob = load_data_ob;
    plot_tracer_ob(obj,out,out_ob,savedir);

    scatter_tracer_ob(obj,out,out_ob,savedir);

end

function out = load_data

    paper_names = {'Wang2020',...
        'Hasenour2020(lac)',...
        'Hasenour2020(prop)',...
        'Hasenour2015',...
        'Satapati2012',...% WT 8wk
        'Burgess2005'};

    flux_names = {'Glucose production',...
        'Glycogenolysis',...
        'GNG from glycerol',...
        'GNG from PEP',...
        'Pyruvate cycling (Pklr)',...
        'Pck1',...
        'TCA cycle'};

    % data_hgp = [60, 110, 109, 60, 2.3/29.3*1000, 0.47*(17.2*4.5/100)/17.2*1000 120];
    data_hgp = [60, 110, 109, 60, 2.3/29.3*1000, 120];

    data_mu = [229.8, 62.9, 114.5*2, 52.4*2, 0, 0, 0;...
        110, 15, 0, 191, 172, 363, 183;...
        109, 14, 0, 189, 82, 272, 195;...
        1, 0.005, 0.256*2, 0.739*2, 1.44/(1.59/0.739/2), 3.03/(1.59/0.739/2), 1/(1.59/0.739/2);...
    %     1, 0.02, 0.28*2, 0.70*2, 2.8/(1.7/0.70/2), 4.5/(1.7/0.70/2), 0/(1.7/0.70/2);...
        1, 0.02, 0.28*2, 0.70*2,[5.31, 8.54, 1.95]./(0.02/2+0.68+1.61);...
        1, 0.11, 0.26*2, 0.63*2, 0.94/(2.12/0.63/2), 3.06/(2.12/0.63/2), 1/(2.12/0.63/2)];

    ci_wang2020 = [sqrt(sum(([68.9 117.5*2 55.4*2]-data_mu(1,2:4)).^2))+data_mu(1,1), 68.9, 117.5*2, 55.4*2, 0, 0, 0]...
        - [-sqrt(sum((-[0.0 105.5*2 43.4*2]+data_mu(1,2:4)).^2))+data_mu(1,1),0.0, 105.5*2 43.4*2 0 0 0];
    ci_wang2020 = ci_wang2020./1.96;

    data_ci_upper = [ci_wang2020;...
        [12 7 0 20 42 59 28];...
        [12 4 0 20 63 74 42];...
        [sqrt(sum(0.004^2+(0.008*2)^2+(0.009*2)^2)), 0.004, 0.008*2, 0.009*2, [0.008, 0.14, 0.09]/(1.59/0.739/2)];...
        [sqrt(sum(0.02^2+(0.02*2)^2+(0.02*2)^2)), 0.02, 0.02*2, 0.02*2, [0.55,0.63,0.15]/((0.02/2+0.68+1.61)/(0.04/2+0.05+0.06)/2)];...
        [sqrt(sum(0.06^2+(0.09*2)^2+(0.06*2)^2)), 0.06, 0.09*2, 0.06*2, [0.46,1.1,0.83]/(2.12/0.63/2)];...
        ];

    data_ci_lower = data_ci_upper;

    out.paper_names = paper_names;
    out.flux_names = flux_names;
    out.data_mu = data_mu./data_mu(:,1);
    out.data_ci_upper = data_ci_upper./data_mu(:,1);
    out.data_ci_lower = data_ci_lower./data_mu(:,1);
    out.data_hgp = data_hgp;

end

function scatter_tracer_ob(obj,out_ctl,out,savedir)
    
    idx_flux = [1,1,4,7,8,12];
    num_flux = length(out.flux_names);
    iter = size(obj.par.a,1);
    v = reshape(obj.par.v(:,1,idx_flux),iter,num_flux);
    v(:,1) = sum(obj.par.v(:,1,[1 2]),[2 3]);
    v(:,3) = sum(obj.par.v(:,1,[4 6]),[2 3]);
    v_ob = reshape(obj.par.v(:,3,idx_flux),iter,num_flux);
    v_ob(:,1) = sum(obj.par.v(:,3,[1 2]),[2 3]);
    v_ob(:,3) = sum(obj.par.v(:,3,[4 6]),[2 3]);
    num_p = length(out.paper_names);

    v_fc = v_ob./v;
    mu_v = median(v_fc,1);
    ci_v = nan(size(mu_v));
    ci_v(1,:) = prctile(v_fc,2.5,1);
    ci_v(2,:) = prctile(v_fc,97.5,1);
%     std_v = std(v_fc,[],1);
%     ci_v = 1.96.*std_v;
%     ci_v = [std(v_fc,[],1);std(v_fc,[],1)];

    v_tracer = out.data_mu;
    col_list = [0 90 255;...% blue
        3 175 122;...% green
        77 196 255;...% skyblue
        ]./255;
    
    ci_data = out.data_ci_lower;
    
    % all studies
    fig = figure('visible','off');
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [4.5 4.5];
    fig.PaperPosition = [0 0 4.5 4.5];
    hold on;
    line([0 4],[0 4],...
        'LineWidth',1,'Color','k')
%     v_tracer_nan = nan(size(v_tracer));
    for p=1:num_p
        v_tracer_tmp = v_tracer(p,:);
%         v_tracer_tmp(v_tracer_tmp==0) = nan;
%         v_tracer_nan(p,:) = v_tracer_tmp;
        eb(1) = errorbar(v_tracer_tmp,mu_v,...
            ci_data(p,:),ci_data(p,:),...
            'horizontal', 'LineStyle', 'none',...
            'Color',col_list(p,:),'LineWidth',1);
        eb(2) = errorbar(v_tracer_tmp,mu_v,...
            mu_v-ci_v(1,:),ci_v(2,:)-mu_v, 'vertical', 'LineStyle', 'none',...
            'Color',col_list(p,:),'LineWidth',1);
        scatter(v_tracer_tmp,mu_v,...
            'SizeData',40,...
            'MarkerFaceColor',col_list(p,:),...
            'MarkerEdgeColor',col_list(p,:));
    end
    v_tmp = repmat(mu_v,num_p,1);
    [cor_p_, pval_p_] = corrcoef(v_tracer(:),v_tmp(:),'Rows','complete');
    [cor_p,pval_p] = corr(v_tracer(:),v_tmp(:),'Type','Pearson','Rows','complete');
    [cor_s, pval_s] = corr(v_tracer(:),v_tmp(:),'Type','Spearman','Rows','complete');
    title(['Pearson = ' num2str(round(cor_p,3)) ' (p='...
        num2str(sprintf('%.2e',pval_p)) ')',...
        ' Spearman = ' num2str(round(cor_s,3)) ' (p='...
        num2str(sprintf('%.2e',pval_s)) ')',...
        ],'FontSize',10);
    xlabel('Fold changes in the previous studies');
    ylabel('Fold changes inferred by OMELET');
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
    print(fig,'-painters',...
        [savedir '/scatter_fluxFC_allGroup.pdf'],...
        '-dpdf','-bestfit');
    close all;

end

function scatter_tracer(obj,out,savedir)
    
    idx_flux = [1,1,4,6,7,8,12];
    num_flux = length(out.flux_names);
    iter = size(obj.par.a,1);
    v = reshape(obj.par.v(:,1,idx_flux),iter,num_flux);
    v(:,1) = sum(obj.par.v(:,1,[1 2]),[2 3]);
    num_p = length(out.paper_names);

    ci_data = out.data_ci_lower;
    
    mu_v = median(v,1);
%     sd_v = std(v,1);
%     ci_v = [sd_v;sd_v];
    mu_v = mu_v./mu_v(1);
%     ci_v = ci_v ./mu_v(1);
    ci_v = nan(2,num_flux);
    ci_v(1,:) = prctile(v,2.5,1);
    ci_v(2,:) = prctile(v,97.5,1);
%     ci_v = [std(v,1);std(v,1)];

    v_tracer = out.data_mu;
    col_list = [150 150 150;...% grey
        0 90 255;...% blue
        77 196 255;...% skyblue
        255 0 0;...% red
        3 175 122;...% green
        246 170 0;...% oragne
        153 0 153;...% purple
        ]./255;
    
    % all studies
    fig = figure('visible','off');
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [4.5 4.5];
    fig.PaperPosition = [0 0 4.5 4.5];
    hold on;
    line([0 4],[0 4],...
        'LineWidth',1,'Color','k')
    v_tracer_nan = nan(size(v_tracer));
    for p=1:num_p
        v_tracer_tmp = v_tracer(p,:);
        v_tracer_tmp(v_tracer_tmp==0) = nan;
        v_tracer_nan(p,:) = v_tracer_tmp;
        eb(1) = errorbar(v_tracer_tmp,mu_v,...
            ci_data(p,:), 'horizontal', 'LineStyle', 'none',...
            'Color',col_list(p,:),'LineWidth',1);
        eb(2) = errorbar(v_tracer_tmp,mu_v,...
            mu_v-ci_v(1,:),ci_v(2,:)-mu_v, 'vertical', 'LineStyle', 'none',...
            'Color',col_list(p,:),'LineWidth',1);
        scatter(v_tracer_tmp,mu_v,...
            'SizeData',40,...
            'MarkerFaceColor',col_list(p,:),...
            'MarkerEdgeColor',col_list(p,:));
    end
    v_tmp = repmat(mu_v,num_p,1);
    [cor_p_, pval_p_] = corrcoef(v_tracer_nan(:),v_tmp(:),'Rows','complete');
    [cor_p,pval_p] = corr(v_tracer_nan(:),v_tmp(:),'Type','Pearson','Rows','complete');
    [cor_s, pval_s] = corr(v_tracer_nan(:),v_tmp(:),'Type','Spearman','Rows','complete');
    title(['Pearson = ' num2str(round(cor_p,3)) ' (p='...
        num2str(sprintf('%.2e',pval_p)) ')',...
        ' Spearman = ' num2str(round(cor_s,3)) ' (p='...
        num2str(sprintf('%.2e',pval_s)) ')',...
        ],'FontSize',10);
    xlabel('Metabolic fluxes in the previous studies');
    ylabel('Metabolic fluxes inferred by OMELET');
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
    print(fig,'-painters',...
        [savedir '/scatter_flux_allGroup.pdf'],...
        '-dpdf','-bestfit');
    close all;
    
    % individual studies
    fig = figure('visible','off');
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [8 5];
    fig.PaperPosition = [0 0 8 5];
    hold on;
    for p=1:num_p
        subplot(2,3,p);
        line([0 4],[0 4],...
            'LineWidth',1,'Color','k')
        hold on;
        v_tracer_tmp = v_tracer(p,:);
        v_tracer_tmp(v_tracer_tmp==0) = nan;
        eb(1) = errorbar(v_tracer_tmp,mu_v,...
            ci_data(p,:), 'horizontal', 'LineStyle', 'none',...
            'Color',col_list(p,:),'LineWidth',1);
        eb(2) = errorbar(v_tracer_tmp,mu_v,...
            mu_v-ci_v(1,:),ci_v(2,:)-mu_v, 'vertical', 'LineStyle', 'none',...
            'Color',col_list(p,:),'LineWidth',1);
        scatter(v_tracer_tmp,mu_v,...
            'SizeData',40,...
            'MarkerFaceColor',col_list(p,:),...
            'MarkerEdgeColor',col_list(p,:));
        text(v_tracer_tmp,mu_v,out.flux_names,'FontSize',10);
        cor_now = corr(v_tracer_tmp',mu_v','rows','complete');
        title(['Pearson = ' num2str(round(cor_now,3))],'FontSize',10);
        xlabel(out.paper_names(p));
        ylabel('OMELET');
        set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
    end
    print(fig,'-painters',...
        [savedir '/scatter_flux_eachGroup.pdf'],...
        '-dpdf','-bestfit');
    close all;

end

function plot_tracer(obj,out,savedir)

    idx_flux = [1,1,4,6,7,8,12];
    num_flux = length(out.flux_names);
    iter = size(obj.par.a,1);
    v = reshape(obj.par.v(:,1,idx_flux),iter,num_flux);
    v(:,1) = sum(obj.par.v(:,1,[1 2]),[2 3]);
    num_p = length(out.paper_names);

    mu_v = median(v,1);
    % sd_v = std(v,1);
    % ci_v = [sd_v;sd_v];
    mu_v = mu_v./mu_v(1);

    ci_v = nan(2,num_flux);
    ci_v(1,:) = prctile(v,2.5,1);
    ci_v(2,:) = prctile(v,97.5,1);
    % ci_v = ci_v ./mu_v(1);


    v_tracer = out.data_mu;
    col_list = [150 150 150;...% grey
        0 90 255;...% blue
        77 196 255;...% skyblue
        255 0 0;...% red
        3 175 122;...% green
        246 170 0;...% oragne
        153 0 153;...% purple
        ]./255;

    bar_input = [mu_v;v_tracer];
    bar_col = [0 0 0; col_list];
    err_input_pos = [ci_v(2,:)-mu_v; out.data_ci_upper];
    err_input_neg = [mu_v-ci_v(1,:); out.data_ci_lower];
    bar_names = out.flux_names;

    %%%%%%%%%%%%%%%%%%
    bar_width = 0.8;
    err_width = 2.0;

    % init defaults for parameters
    [N_grps,N_bars]=size(bar_input);

    % init group width and bar shift
    shift_span=(1-bar_width)*(N_grps-1);
    bar_shift=linspace(-shift_span/2,+shift_span/2,N_grps);

    % init handles vectors
    hb=zeros(N_grps,1);
    he=zeros(N_grps,1);

    fig = figure('visible','off');
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    hold on;
    % plot bars
    for grp=1:N_grps
        xpos = (grp:N_grps:N_bars*N_grps-(N_grps-grp))-bar_shift(grp);
        hb(grp) = bar(xpos,bar_input(grp,:),...
            'BarWidth',bar_width/N_grps,...
            'FaceColor',bar_col(grp,:),...
            'EdgeColor','none'); 
    end

    for grp=1:N_grps
        xpos = (grp:N_grps:N_bars*N_grps-(N_grps-grp))-bar_shift(grp);
        err_lower = err_input_neg(grp,:);
        err_upper = err_input_pos(grp,:);
       he(grp) = errorbar(xpos,bar_input(grp,:),err_lower,err_upper,...
           'CapSize',err_width,...
           'Color','k',...
           'Marker','none',...
           'LineStyle','none'); 
    end

    legend(['OMELET', out.paper_names],'Location','southeastoutside');

    % compute position of group x ticks
    bar_xtick=N_grps/2+0.5:N_grps:N_bars*N_grps-N_grps/2+0.5;
    % set the x tick labels
    set(gca,'XTick',bar_xtick,'XTickLabel',bar_names);
    % cosmetic fine-tuning of the figure
    set(gca,'XLim',[0 bar_xtick(end)+bar_xtick(1)]); % adjusts the x axis to the plot

    set(gca,'XTickLabelRotation',30);
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
    ylabel('Relative flux value');

    fig.PaperSize = [8 3];
    fig.PaperPosition = [0 0 8 3];
    % saveas( gcf, [ savedir '/v_tracer_cmp_WT0h.png' ] );
    saveas( gcf, [ savedir '/v_tracer_cmp_WT0h.pdf' ] );
    close all;

    % output as table
    stat_names = {'Mean','SD'};
    colnames = cellfun(@(y) cellfun(@(x) [y ' (' x ')'], stat_names,'UniformOutput',false),...
        ['OMELET' out.paper_names],'UniformOutput',false);
    tmp = {};
    for i=1:num_p+1
        tmp = [tmp; colnames{i}; num2cell([bar_input(i,:)' err_input_pos(i,:)'])];
    end
    rownames = repmat([{''} out.flux_names]',num_p+1,1);

    tbl_out = [[{''};rownames] [stat_names;tmp]];
    % tbl_out = [tbl_out ['Absolute glucose production (umol/kg/min)'; num2cell([mean(out.data_hgp) out.data_hgp])']];
    writecell(tbl_out,[savedir '/data_cmpTracer.csv']);

end

function out = load_data_ob

    paper_names = {'Patterson2016',...
        'Satapati2012',...
        'Turner2005'};

    flux_names = {'Glucose production',...
        'Glycogenolysis',...
        'GNG',...
        'Pyruvate cycling (Pklr)',...
        'Pck1',...
        'TCA cycle'};

    data_mu_ctl = [4.2 0 [0 4.9 9.1 2.6].*2;...
        2.4 0.14 [(0.47+2.28) 4.65 8.29 2.71].*2;...
        25.8 4.1-0.3 [21.7*2 0 0 0].*2];
    data_mu_ob = [5.9 0 [0 11.1 16.9 5.4].*2;...
        3.6 0.51 [0.64+3.03 8.27 13.1 3.71].*2;...
        44.7 12.5-3.2 32.1*2 0 0 0];
    data_mu = data_mu_ob./data_mu_ctl;

    % coef_ctl = [1000/30.9 1000/38.5 1000*(10.6/12.1)/180];
    % coef_ob = [1000/37.2 1000/53.6 1000*(10.7/24.8)/180];

    data_hgp_ctl = [4.2/30.9*1000 2.4/38.5*1000, 25.8*(10.6/12.1)/180*1000 ];
    data_hgp_ob = [5.9/37.2*1000 3.6/53.6*1000, 44.7*(10.7/24.8)/180*1000 ];


    % data_mu = (data_mu_ob.*coef_ob')./(data_mu_ctl.*coef_ctl');

    % CI = 1.96.*data_sd./sqrt(data_n) = 1.96.*data_se;
    data_ctl_ci_upper = [0.4 0 0 1.0 1.2 0.5;...
        0.1 0.05 0.19 0.34 0.31 0.28;...
        6.0 1.0 5.8 0 0 0];
    data_ob_ci_upper = [0.4 0 0 1.9 2.2 0.6;...
        0.3 0.07 0.32 0.95 1.15 0.45;...
        11.3 3.8 7.1 0 0 0];
    data_ci_upper = sqrt((data_ob_ci_upper./data_mu_ctl).^2 + ...
        (data_ctl_ci_upper.*data_mu_ob./(data_mu_ctl).^2).^2);

    out.paper_names = paper_names;
    out.flux_names = flux_names;
    out.data_mu = data_mu;
    out.data_ci_upper = data_ci_upper;
    out.data_ci_lower = data_ci_upper;
    out.data_ci_upper_ctl = data_ctl_ci_upper;
    out.data_ci_upper_ob = data_ob_ci_upper;
    out.data_ci_lower_ctl = data_ctl_ci_upper;
    out.data_ci_lower_ob = data_ob_ci_upper;
    out.data_mu_ctl = data_mu_ctl;
    out.data_mu_ob = data_mu_ob;
    out.data_hgp_ctl = data_hgp_ctl;
    out.data_hgp_ob = data_hgp_ob;

end

function plot_tracer_ob(obj,out_ctl,out,savedir)

    idx_flux = [1,1,4,7,8,12];
    num_flux = length(out.flux_names);
    iter = size(obj.par.a,1);
    v = reshape(obj.par.v(:,1,idx_flux),iter,num_flux);
    v(:,1) = sum(obj.par.v(:,1,[1 2]),[2 3]);
    v(:,3) = sum(obj.par.v(:,1,[4 6]),[2 3]);
    v_ob = reshape(obj.par.v(:,3,idx_flux),iter,num_flux);
    v_ob(:,1) = sum(obj.par.v(:,3,[1 2]),[2 3]);
    v_ob(:,3) = sum(obj.par.v(:,3,[4 6]),[2 3]);
    num_p = length(out.paper_names);

    v_fc = v_ob./v;
    mu_v = median(v_fc,1);
    % std_v = std(v_fc,[],1);
    % ci_v = std_v;
    ci_v = nan(size(mu_v));
    ci_v(1,:) = prctile(v_fc,2.5,1);
    ci_v(2,:) = prctile(v_fc,97.5,1);

    v_tracer = out.data_mu;
    col_list = [150 150 150;...% grey
        0 90 255;...% blue
        77 196 255;...% skyblue
        255 0 0;...% red
        3 175 122;...% green
        246 170 0;...% oragne
        153 0 153;...% purple
        ]./255;

    bar_input = [mu_v;v_tracer];
    bar_col = [0 0 0; col_list([2 5 3],:)];
    err_input_pos = [ci_v(2,:)-mu_v; out.data_ci_upper];
    err_input_neg = [mu_v-ci_v(1,:); out.data_ci_lower];
    bar_names = out.flux_names;

    %%%%%%%%%%%%%%%%%%
    bar_width = 0.8;
    err_width = 4.0;

    % init defaults for parameters
    [N_grps,N_bars]=size(bar_input);

    % init group width and bar shift
    shift_span=(1-bar_width)*(N_grps-1);
    bar_shift=linspace(-shift_span/2,+shift_span/2,N_grps);

    % init handles vectors
    hb=zeros(N_grps,1);
    he=zeros(N_grps,1);

    fig = figure('visible','off');
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    hold on;
    % plot bars
    for grp=1:N_grps
        xpos = (grp:N_grps:N_bars*N_grps-(N_grps-grp))-bar_shift(grp);
        hb(grp) = bar(xpos,bar_input(grp,:),...
            'BarWidth',bar_width/N_grps,...
            'FaceColor',bar_col(grp,:),...
            'EdgeColor','none'); 
    end

    for grp=1:N_grps
        xpos = (grp:N_grps:N_bars*N_grps-(N_grps-grp))-bar_shift(grp);
        err_lower = err_input_neg(grp,:);
        err_upper = err_input_pos(grp,:);
       he(grp) = errorbar(xpos,bar_input(grp,:),err_lower,err_upper,...
           'CapSize',err_width,...
           'Color','k',...
           'Marker','none',...
           'LineStyle','none'); 
    end

    legend(['OMELET', out.paper_names],'Location','southeastoutside');

    % compute position of group x ticks
    bar_xtick=N_grps/2+0.5:N_grps:N_bars*N_grps-N_grps/2+0.5;
    % set the x tick labels
    set(gca,'XTick',bar_xtick,'XTickLabel',bar_names);
    % cosmetic fine-tuning of the figure
    set(gca,'XLim',[0 bar_xtick(end)+bar_xtick(1)]); % adjusts the x axis to the plot

    set(gca,'XTickLabelRotation',30);
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
    ylabel('Flux fold change in obesity');

    fig.PaperSize = [6 3];
    fig.PaperPosition = [0 0 6 3];
    % saveas( gcf, [ savedir '/v_tracer_cmp_obFC.png' ] );
    saveas( gcf, [ savedir '/v_tracer_cmp_obFC.pdf' ] );
    close all;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v_ob_abs = mean(v_ob,1) .* mean(out_ctl.data_hgp);
    v_ctl_abs = v .* mean(out_ctl.data_hgp);
    v_tracer_ob = (out.data_mu_ob./out.data_mu_ob(:,1)) .* repmat(out.data_hgp_ob',1,num_flux);
    v_tracer_ctl = out.data_mu_ctl;
    bar_input = [v_ob_abs;v_tracer_ob];

    err_input_neg = [std(v_ob,[],1).*1.96; out.data_ci_lower_ob./out.data_mu_ob(:,1)] .* repmat([mean(out_ctl.data_hgp) out.data_hgp_ob]',1,num_flux);
    err_input_pos = [std(v_ob,[],1).*1.96; out.data_ci_upper_ob./out.data_mu_ob(:,1)] .* repmat([mean(out_ctl.data_hgp) out.data_hgp_ob]',1,num_flux);

    bar_width = 0.8;
    err_width = 4.0;

    % init defaults for parameters
    [N_grps,N_bars]=size(bar_input);

    % init group width and bar shift
    shift_span=(1-bar_width)*(N_grps-1);
    bar_shift=linspace(-shift_span/2,+shift_span/2,N_grps);

    % init handles vectors
    hb=zeros(N_grps,1);
    he=zeros(N_grps,1);



    % output as table
    stat_names = {'Mean','SD'};
    colnames = cellfun(@(y) cellfun(@(x) [y ' (' x ')'], stat_names,'UniformOutput',false),...
        ['OMELET' out.paper_names],'UniformOutput',false);
    tmp = {};
    for i=1:num_p+1
        tmp = [tmp; colnames{i}; num2cell([bar_input(i,:)' err_input_pos(i,:)'])];
    end
    rownames = repmat([{''} out.flux_names]',num_p+1,1);

    tbl_out = [[{''};rownames] [stat_names;tmp]];
    % tbl_out = [tbl_out ['Absolute glucose production (umol/kg/min)'; num2cell([mean(out.data_hgp) out.data_hgp])']];
    writecell(tbl_out,[savedir '/data_cmpTracer_ob.csv']);



end
