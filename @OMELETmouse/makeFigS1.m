function makeFigS1(obj,savedir)
    
    savedir = [savedir '/FigS1'];
    mkdir(savedir);
    
    blood = obj.omics_data.blood;

    % plasma glucose, plasma insulin and glycogen
    % plot_plasma(blood,'glucose',savedir);
    % plot_plasma(blood,'insulin',savedir);
    plot_plasma(obj,blood,'glycogen',savedir);

    % time-course of plasma glucose and insulin
    dat_csv = 'plasma_glucose';
    plot_timecourse(obj,dat_csv,savedir);

    dat_csv = 'plasma_insulin';
    plot_timecourse(obj,dat_csv,savedir);
    
    % PCA
    plot_PCA(obj,savedir);
    
end

function plot_plasma(obj,blood,var_str,savedir)

    info_grp = obj.model_data;

    switch var_str
        case 'glycogen'
            data = [blood.glyc0min([13:23 1:12]) blood.glyc240min([13:24 1:12])];
            grp_names = info_grp.grp_names;
            num_grp = info_grp.num_smpl_g;
            idx_grp = info_grp.idx_g;
            col = info_grp.col;

        case 'insulin'
            data = blood.insulin0min([12:23 1:11]);
            grp_names = info_grp.grp_names([1 3]);
            num_grp = info_grp.num_smpl_g([1 3]);
            idx_grp = info_grp.idx_g([1 2],:);
            col = info_grp.col([1 3],:);

        case 'glucose'
            data = blood.bs0min([12:23 1:11]);
            grp_names = info_grp.grp_names([1 3]);
            num_grp = info_grp.num_smpl_g([1 3]);
            idx_grp = info_grp.idx_g([1 2],:);
            col = info_grp.col([1 3],:);

    end


    c = categorical(grp_names,grp_names);

    % plot per group
    fig = figure('visible','off');
    for g=1:length(grp_names)
        scatter(repmat(c(g),num_grp(g),1),data(idx_grp(g,1):idx_grp(g,2)),...
            20,col(g,:),'filled');
        hold on;
    end
    ylabel(var_str);
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 length(grp_names) 3];
    fig.PaperSize = [length(grp_names) 3];
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
%     saveas(gcf,[savedir '/scatter_' var_str '.png' ]);
    saveas(gcf,[savedir '/scatter_' var_str '.pdf' ]);
    close all;

    % plot by mean and sd
    mu = nan(1,length(grp_names));
    sd = nan(1,length(grp_names));
    for g=1:length(grp_names)
        mu(g) = mean(data(idx_grp(g,1):idx_grp(g,2)));
        sd(g) = std(data(idx_grp(g,1):idx_grp(g,2)));
    end
    dat_mu = reshape(mu,2,2)';
    dat_sd = reshape(sd,2,2)';
    col = [info_grp.col(1,:); info_grp.col(3,:)];
    t = [0 240];

    fig = figure('visible','off');
    hold on;
    for i=1:2
        line(t,dat_mu(i,:),'LineWidth',2,'Color',col(i,:));
        errorbar(t,dat_mu(i,:),dat_sd(i,:),'Color',col(i,:));
    end
    ax = gca;
    ax.XTick = t;
    ax.XLim = [t(1)-10 t(2)+10];
    ylabel(var_str);
    xlabel('time(min)');

    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 2 3];
    fig.PaperSize = [2 3];
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
%     saveas(gcf,[savedir '/line_' var_str '.png' ]);
    saveas(gcf,[savedir '/line_' var_str '.pdf' ]);
    close all;

end

function plot_timecourse(obj,dat_csv,savedir)
 
    info_grp = obj.model_data;
    
    x = readtable([obj.data_path '/' dat_csv '.csv']);
    t = x{1,2:end};
    varname = x{2:end,1};
    grp_idx = contains(varname,'WT');
    dat = x{2:end,2:end};
    
    dat_mu = [mean(dat(grp_idx,:),1); mean(dat(~grp_idx,:),1) ];
    dat_std = [std(dat(grp_idx,:),[],1); std(dat(~grp_idx,:),[],1)];
    
    col = [info_grp.col(1,:); info_grp.col(3,:)];
    mouse = {'WT','ob'};
    
    % both
    fig = figure('visible','off');
    hold on;
    for i=1:2
        line(t,repmat(dat_mu(i,1),1,length(t)),...
            'LineWidth',2,'Color',col(i,:),'LineStyle','--');
        line(t,dat_mu(i,:),'LineWidth',2,'Color',col(i,:));
        errorbar(t,dat_mu(i,:),dat_std(i,:));
    end
    xticks([0 60 120 180 240]);
    xlabel('time(min)');
    switch dat_csv
        case 'plasma_glucose'
            ylab_now = 'mg/dl';
        case 'plasma_insulin'
            ylab_now = 'pg/ml';
    end
    ylabel(ylab_now);
    title(dat_csv,'Interpreter','none');
    
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 4 3];
    fig.PaperSize = [4 3];
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
%     saveas(gcf,[savedir '/timecourse_' dat_csv '.png' ]);
    saveas(gcf,[savedir '/timecourse_' dat_csv '.pdf' ]);
    close all;
    
    % separately
    fig = figure('visible','off');
    hold on;
    for i=1:2
        subplot(2,1,i);
        hold on;
        line(t,repmat(dat_mu(i,1),1,length(t)),...
            'LineWidth',2,'Color',col(i,:),'LineStyle','--');
        line(t,dat_mu(i,:),'LineWidth',2,'Color',col(i,:));
        errorbar(t,dat_mu(i,:),dat_std(i,:),'Color',col(i,:));
        xticks([0 60 120 180 240]);
        xlabel('time(min)');
        ylabel(ylab_now);
        title(mouse{i});
    end
    sgtitle(dat_csv,'Interpreter','none');
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 4 6];
    fig.PaperSize = [4 6];
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
%     saveas(gcf,[savedir '/timecourse_subplot_' dat_csv '.png' ]);
    saveas(gcf,[savedir '/timecourse_subplot_' dat_csv '.pdf' ]);
    close all;
    
end

function plot_PCA(obj,savedir)
    
    data_omics = obj.omics_data.data_omics;
    dict_omics = obj.omics_data.dict_omics;
    name_omics = obj.omics_data.name_omics;
    info_grp = obj.model_data;
    
    data_all = [];
    var_all = {};
    idx_var = nan(length(data_omics),2);
    idx_tmp = 1;
    for i=1:length(data_omics)
        data_all = [data_all data_omics{i}];
        var_all = [var_all; dict_omics{i}.val(:)];
        idx_var(i,1) = idx_tmp;
        idx_var(i,2) = idx_tmp + length(dict_omics{i}.val)-1;
        idx_tmp = idx_var(i,2)+1;
    end

    data_all_n = data_all ./ nanmean(data_all,1);
    num_var = size(data_all_n,2);

    % remove nan vars
    idx_nan = false(info_grp.num_g,num_var);
    for g=1:info_grp.num_g
        data_now = data_all_n(info_grp.idx_g(g,1):info_grp.idx_g(g,2),:);
        idx_nan(g,:) = any(isnan(data_now),1);
    end
    idx_analyze = ~any(idx_nan,1);
    data_input = data_all_n(:,idx_analyze);
    var_input = var_all(idx_analyze);

    idx_var_input = nan(length(data_omics),2);
    idx_tmp = 1;
    for i=1:length(data_omics)
        idx_var_input(i,1) = idx_tmp;
        num_tmp = sum(idx_analyze(idx_var(i,1):idx_var(i,2)));
        idx_var_input(i,2) = idx_tmp + num_tmp-1;
        idx_tmp = idx_var_input(i,2)+1;
    end
    assert(idx_var_input(end,end)==length(var_input));

    % PCA for each omics
    for i=1:length(data_omics)
        [coeff,score,latent,tsquared,explained,mu] =...
            pca(data_input(:,idx_var_input(i,1):idx_var_input(i,2)));
        % score plot
        plot_score(score,explained,info_grp,name_omics{i},savedir);
    end


    [coeff,score,latent,tsquared,explained,mu] = pca(data_input);

    % score plot
    plot_score(score,explained,info_grp,'all',savedir);

    % loading plot
    % plot_loading(coeff,var_input,idx_var_input,savedir);

    % plot_score_loading(coeff,score,explained,info_grp,var_input,idx_var_input,savedir)

%     save([savedir '/workspace_figS2_' datestr(now,'yyyymmdd') '.mat']);

    
end


function plot_score(score,explained,info,types,savedir)

% PC pairs to plot
cmb_pc = [1 2;1 3;2 3];

for c=1:size(cmb_pc,1)
    pc_1 = cmb_pc(c,1);
    pc_2 = cmb_pc(c,2);
    fig = figure('visible','off');
    hold on;
    for g=1:info.num_g
        idx_now = info.idx_g(g,1):info.idx_g(g,2);
        scatter(score(idx_now,pc_1),score(idx_now,pc_2),...
            20,'filled','CData',info.col(g,:));
    end
    xlabel(['PC' num2str(pc_1) ' (' num2str(round(explained(pc_1),2)) '%)']);
    ylabel(['PC' num2str(pc_2) ' (' num2str(round(explained(pc_2),2)) '%)']);
    legend(info.grp_names,'Location','southeastoutside');

    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [ 5 4 ];
    fig.PaperPosition = [ 0 0 5 4 ];
    set( findobj(gcf, 'Type','Axes'), 'FontSize', 10 );
%     saveas(gcf,[savedir '/pca_score_',types,'_pc' num2str(pc_1) '_pc' num2str(pc_2) '.png']);
    saveas(gcf,[savedir '/pca_score_',types,'_pc' num2str(pc_1) '_pc' num2str(pc_2) '.pdf']);
    close all;
end

end