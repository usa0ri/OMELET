function makeFig2(obj,savedir_)
    
    savedir = [savedir_ '/Fig2'];
    mkdir(savedir);
    
   obj.model_data.col = obj.model_data.col;
   var_omics = obj.omics_data.var_omics;
   data_omics = obj.omics_data.data_omics;
   name_omics = obj.omics_data.name_omics;
                        
    pmat = cell(1,3);
    for i=1:length(var_omics)
        plot_indv_bar_omics(data_omics{i},var_omics{i},name_omics{i},obj.model_data,savedir);
        plot_bar_omics(data_omics{i},var_omics{i},name_omics{i},obj.model_data,savedir);
        pmat{i} = ttest_omics(data_omics{i},var_omics{i},name_omics{i},obj.model_data,savedir);
    end
    obj.omics_data.pmat = pmat;
    
end


function plot_bar_omics(data_now,var_now,name_now,info,savedir)

    [num_smpl, num_var] = size(data_now);

    % fig = figure('visible','off');
    mu = nanmean(data_now(info.idx_g(1,1):info.idx_g(1,2),:),1);
    data_n = data_now./ mu;
    bar_input = nan(info.num_g,num_var);
    err_input = nan(info.num_g,num_var);
    for g=1:info.num_g
        bar_input(g,:) = nanmean(data_n(info.idx_g(g,1):info.idx_g(g,2),:),1);
        err_input(g,:) = nanstd(data_n(info.idx_g(g,1):info.idx_g(g,2),:),[],1);
    end
    % assert(all(bar_input(1,:)-1<1e-3));

    bar_col = info.col;
    bar_names = var_now;
    
    %%%%%%%%%%%%%%
    % reorder sample index
    idx_s = info.idx_g_reorder;
    bar_input = bar_input(idx_s,:);
    err_input = err_input(idx_s,:);
    bar_col = bar_col(idx_s,:);
    %%%%%%%%%%%%%%
    

    plot_errbars(bar_input,err_input,bar_col,bar_names,name_now,savedir);

end

function plot_indv_bar_omics(data_now,var_now,name_now,info,savedir)

    [num_smpl, num_var] = size(data_now);

    % fig = figure('visible','off');
    mu = nanmean(data_now(info.idx_g(1,1):info.idx_g(1,2),:),1);
    data_n = data_now./ mu;
    bar_input = nan(info.num_g,num_var);
    errorbar_input = nan(info.num_g,num_var);
    for g=1:info.num_g
        bar_input(g,:) = nanmean(data_n(info.idx_g(g,1):info.idx_g(g,2),:),1);
        errorbar_input(g,:) = nanstd(data_n(info.idx_g(g,1):info.idx_g(g,2),:),[],1);
    end
    % assert(all(bar_input(1,:)-1<1e-3));

    %%%%%%%%%%%%%%
    % reorder sample index
    idx_s = info.idx_g_reorder;
    bar_input = bar_input(idx_s,:);
    errorbar_input = errorbar_input(idx_s,:);
    bar_col = info.col(idx_s,:);
    %%%%%%%%%%%%%%
    
    c = categorical(info.grp_names,info.grp_names);

    fig = figure('visible','off');
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    num_sub = ceil(sqrt(num_var));
    for i=1:num_var
        subplot(num_sub,num_sub,i);
        b = bar(c,bar_input(:,i),'FaceColor','flat');
        b.CData = bar_col; 
        b.BarWidth = 0.8;
        b.EdgeColor = 'none';
        hold on;
        e = errorbar(c,bar_input(:,i),errorbar_input(:,i),'.k');
        e.LineWidth = 0.5;
        e.Marker = 'none';
        title(var_now{i});
        ax = gca;
        ax.XTickLabel = '';
        set(findobj(gca,'type','axes'),'FontSize',8,'FontName','SansSerif');
    end
    fig.PaperPosition = [0 0 num_sub num_sub];
    fig.PaperSize = [num_sub num_sub];
    % saveas(gcf,[savedir '/plot_bar_indv_' name_now '.png' ]);
    saveas(gcf,[savedir '/plot_bar_indv_' name_now '.pdf' ]);

    close all;

end

function plot_errbars(bar_input,err_input,bar_col,bar_names,name_now,savedir)

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
    if ~isempty(err_input)
        for grp=1:N_grps
            xpos = (grp:N_grps:N_bars*N_grps-(N_grps-grp))-bar_shift(grp);
            err_lower = err_input(grp,:);
            err_upper = err_input(grp,:);
           he(grp) = errorbar(xpos,bar_input(grp,:),err_lower,err_upper,...
               'CapSize',err_width,...
               'Color','k',...
               'Marker','none',...
               'LineStyle','none'); 
        end
    end

    % compute position of group x ticks
    bar_xtick=N_grps/2+0.5:N_grps:N_bars*N_grps-N_grps/2+0.5;
    % set the x tick labels
    set(gca,'XTick',bar_xtick,'XTickLabel',bar_names);
    % cosmetic fine-tuning of the figure
    set(gca,'XLim',[0 bar_xtick(end)+bar_xtick(1)]); % adjusts the x axis to the plot

    sz1 = N_bars/3;
    sz2 = N_bars/6;
    fig.PaperPosition = [0 0 sz1 sz2];
    fig.PaperSize = [sz1 sz2];

    ylabel(['Normalized ' name_now ' amount']);
    set(gca,'XTickLabelRotation',90);
    set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
    % saveas(gcf,[savedir '/plot_bar_' name_now '.png' ]);
    saveas(gcf,[savedir '/plot_bar_' name_now '.pdf' ]);

    close all;

end

function p_mat = ttest_omics(data_now,var_now,name_now,info,savedir)

    var_names = var_now;
    grp_names = info.grp_names;
    idx_cmp = [1 3;2 4;1 2;3 4];

    p_mat = nan(size(idx_cmp,1)*2,length(var_names));
    for c=1:size(idx_cmp,1)
        grp1 = grp_names{idx_cmp(c,1)};
        grp2 = grp_names{idx_cmp(c,2)};
        idx_c1 = info.idx_g(idx_cmp(c,1),:);
        idx_c2 = info.idx_g(idx_cmp(c,2),:);
        data_c1 = data_now(idx_c1(1):idx_c1(2),:);
        data_c2 = data_now(idx_c2(1):idx_c2(2),:);
        p_list = [];
        for i=1:length(var_names)
            [h,p,ci,stats] = ttest2(data_c1(:,i),data_c2(:,i),'VarType','unequal');
            p_list = [p_list p];
        end
        q = mafdr(p_list,'BHFDR',true);
        p_list = q;
        is_sig1 = p_list<0.05;
    %     is_sig2 = p_list<0.01;

        % output
        p_mat(2*(c-1)+1,:) = p_list;
        p_mat(2*c,:) = is_sig1;

        % volcano plot
        fc = nanmean(data_c2,1)./nanmean(data_c1,1);
        thres_fc = 1.5;
        thres_p = 0.05;
        fig = figure('visible','off');
        hold on;
        line(repmat(log2(thres_fc),1,2),[floor(min(-log10(p_list))) ceil(max(-log10(p_list)))],...
            'LineStyle','--','Color','b');
        line(repmat(-log2(thres_fc),1,2),[floor(min(-log10(p_list))) ceil(max(-log10(p_list)))],...
            'LineStyle','--','Color','b');
        line([floor(min(log2(fc))) ceil(max(log2(fc)))], repmat(-log10(thres_p),1,2),...
            'LineStyle','--','Color','b');
        scatter(log2(fc),-log10(p_list),'filled','SizeData',10,'CData',[0.2 0.2 0.2]);
        text(log2(fc)+0.01,-log10(p_list),var_names,...
            'FontSize',10);
        title([name_now '(' grp2 '/' grp1 ')'])

        fig.PaperPosition = [0 0 8 8];
        fig.PaperSize = [8 8];

        ylabel(['-log10(FDR-adjusted q value)']);
        xlabel(['log2(fold change) (' grp2 '/' grp1 ')' ]);
        set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
    %     saveas(gcf,[savedir '/volcano_' name_now '_' grp1 grp2 '.png' ]);
        saveas(gcf,[savedir '/volcano_' name_now '_' grp1 grp2 '.pdf' ]);

        close all;

    end

    % output
    grp_names_cmp = cell(1,size(idx_cmp,1)*2);
    for c=1:size(idx_cmp,1)
        grp_names_cmp{2*(c-1)+1} = [grp_names{idx_cmp(c,1)} grp_names{idx_cmp(c,2)} ];
        grp_names_cmp{2*c} = 'p<0.05';
    end

    out = [var_names num2cell(p_mat')];
    out = [[{''} grp_names_cmp]; out];
    T = cell2table(out);
    writetable(T,[savedir '/pmat_' name_now '.csv'],...
        'WriteRowNames',false,'WriteVariableNames',false);

end
