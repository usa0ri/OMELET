function calcDiffFlux(obj,savedir)

    mkdir(savedir);

    model_data = obj.model_data;
    par = obj.par;

    fluxnames = model_data.X.rxn.rxn_names_include;
    num_rc = model_data.X.num.num_rc;
    iter = size(par.a,1);
    grp_names = model_data.grp_names;
    num_g = length(model_data.grp_names);

    %%%%%%%%%%%%%%
    % reorder fluxnames
    idx_ = obj.fig_info.idx_reorder;
    % fluxnames = fluxnames(idx_);
    rxn_names = obj.fig_info.fluxnames;
    v_tmp = par.v(:,:,idx_);
    %%%%%%%%%%%%%%
    
    % except reactions in TCA cycle
    rxn_names = rxn_names(1:num_rc-5);
    v_tmp = v_tmp(:,:,1:num_rc-5);
    num_rc = num_rc-5;
    
    cmb_list = nchoosek(1:num_rc,2);
    num_c = size(cmb_list,1);
    v_diff = nan(size(v_tmp,1),num_g,size(cmb_list,1));
    rowname_tmp = cell(size(cmb_list,1),1);
    for i=1:size(cmb_list,1)
        rowname_tmp{i} = [rxn_names{cmb_list(i,2)} ' - ' rxn_names{cmb_list(i,1)}];
        for ii=1:num_g
            v_diff(:,ii,i) = v_tmp(:,ii,cmb_list(i,2))-v_tmp(:,ii,cmb_list(i,1));
        end
    end

    tmp = 0.3;
    for g=1:num_g
        v_now = reshape(v_diff(:,g,:),iter,num_c);

        fig  = figure('visible','off');
        boxplot(v_now,...
            'labels',rowname_tmp,...
            'Colors',[0.2 0.2 0.2],...
            'symbol','',...
            'OutlierSize',3,...
            'Width',0.5);

        % get the necessary handles.
        b = get(gca,'children');
        if length(b)>1
            bb = get(b(end),'children'); 
        else
            bb = get(b,'children'); 
        end
        UW = findobj(bb,'Tag','Upper Whisker');
        UAJ = findobj(bb,'Tag','Upper Adjacent Value');
        OUT = findobj(bb,'Tag','Outliers');
        LW = findobj(bb,'Tag','Lower Whisker');
        LAJ = findobj(bb,'Tag','Lower Adjacent Value');

        percen = 97.5;
        for i=1:num_c
            data_sub = v_now(:,num_c-i+1);
            data_highpercent(i) = prctile(data_sub,percen,1);
            data_lowpercent(i) = prctile(data_sub,100-percen,1);
            UW_lims = get(UW(i),'ydata');
            set(UW(i),'ydata',[UW_lims(1) data_highpercent(i)])
            set(UAJ(i),'ydata',[data_highpercent(i) data_highpercent(i)])
            LW_lims = get(LW(i),'ydata');
            set(LW(i),'ydata',[data_lowpercent(i) LW_lims(2)])
            set(LAJ(i),'ydata',[data_lowpercent(i) data_lowpercent(i)])
            outliers = get(OUT(i),'ydata');
            outlier_index = find(outliers == data_highpercent(i));
            outliers(1:outlier_index) = data_highpercent(i);
            set(OUT(i),'ydata',outliers);      
        end
        
        if g<3
            ylim([-1.5 3]);
        else
            ylim([-2.5 5.5]);
        end
        ax = gca;
        ax.XTickLabelRotation = 90;
        hold on;
        line(linspace(ax.XLim(1),ax.XLim(2),10),zeros(1,10),...
            'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
        ylabel(['Flux difference (' grp_names{g} ')']);
        set(findobj(gca,'type','line'),'linew',1);
        set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
        
        idx_inc = data_highpercent>0 & data_lowpercent>0;
        idx_dec = data_highpercent<0 & data_lowpercent<0;
        str_color = repmat({'black'},num_c,1);
        str_color(fliplr(idx_inc)) = {'magenta'};
        str_color(fliplr(idx_dec)) = {'cyan'};
        label_col = strcat('\color{',str_color,'}',rowname_tmp);
        ax.XTickLabel = label_col;
        ax.TickLabelInterpreter = 'tex';
%         ax.TickLength = [0 0];

        fig.PaperUnits = 'inches';
        fig.PaperSize = [10 8];
        fig.PaperPosition = [0 0 10 8];
    %     saveas( gcf, [ savedir '/v_fc_' grp_names{cmb(c,2)} '_' grp_names{cmb(c,1)} '.png' ] );
        saveas( gcf, [ savedir '/v_diff_' grp_names{g} '.pdf' ] );
        close all;

    end

end
