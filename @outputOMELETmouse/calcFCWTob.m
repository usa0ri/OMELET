function calcFCWTob(obj,savedir)

    mkdir(savedir);

    model_data = obj.model_data;
    par = obj.par;

    fluxnames = model_data.X.rxn.rxn_names_include;
    num_rc = model_data.X.num.num_include;
    iter = size(par.a,1);

    %%%%%%%%%%%%%%
    % reorder sample indx
    idx_s = model_data.idx_g_reorder;
    cmb = obj.cont.cmb;
    num_c = size(cmb,1);
    grp_names = model_data.grp_names;

    % reorder fluxnames
    idx_ = obj.fig_info.idx_reorder;
    % fluxnames = fluxnames(idx_);
    fluxnames = obj.fig_info.fluxnames;
    par.v = par.v(:,idx_s,idx_);
    %%%%%%%%%%%%%%

    tmp = 0.3;
    fc = nan(iter,num_rc,num_c);
    for c=1:num_c
        fc_now = reshape(par.v(:,cmb(c,2),:)./par.v(:,cmb(c,1),:) ,iter,num_rc);

        fig  = figure('visible','off');
        boxplot(fc_now,...
            'labels',fluxnames,...
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
        for i=1:num_rc
            data_sub = fc_now(:,num_rc-i+1);
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

        ylim([0 4]);
        ax = gca;
        ax.XTickLabelRotation = 90;
        hold on;
        line(linspace(ax.XLim(1),ax.XLim(2),10),ones(1,10),...
            'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
        ylabel(['FC (' grp_names{cmb(c,2)} '/' grp_names{cmb(c,1)} ')']);
        set(findobj(gca,'type','line'),'linew',1);
        set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');

        fig.PaperUnits = 'inches';
        fig.PaperSize = [num_rc*tmp num_rc*tmp];
        fig.PaperPosition = [0 0 num_rc*tmp num_rc*tmp];
    %     saveas( gcf, [ savedir '/v_fc_' grp_names{cmb(c,2)} '_' grp_names{cmb(c,1)} '.png' ] );
        saveas( gcf, [ savedir '/v_fc_' grp_names{cmb(c,2)} '_' grp_names{cmb(c,1)} '.pdf' ] );
        close all;

        fc(:,:,c) = fc_now;
    end

    obj.fc_wtob = fc;
    idx_cmb = model_data.idx_cmb;
    col_list = model_data.col_cmb;

    for c=1:size(idx_cmb,1)

        fig  = figure('visible','off');
        hold on;
        for ii=1:size(idx_cmb,2)
            fc_now = reshape(par.v(:,cmb(idx_cmb(c,ii),2),:)./par.v(:,cmb(idx_cmb(c,ii),1),:) ,iter,num_rc);
            boxplot(fc_now,...
                'labels',fluxnames,...
                'Colors',col_list(2*(c-1)+ii,:),...
                'symbol','',...
                'OutlierSize',3,...
                'Width',0.25,...
                'Positions',(1:length(fluxnames))+0.125*(-1)^ii);

            % get the necessary handles.
            b = get(gca,'children');
            if length(b)>1
                bb = get(b(1),'children'); 
            else
                bb = get(b,'children'); 
            end
            UW = findobj(bb,'Tag','Upper Whisker');
            UAJ = findobj(bb,'Tag','Upper Adjacent Value');
            OUT = findobj(bb,'Tag','Outliers');
            LW = findobj(bb,'Tag','Lower Whisker');
            LAJ = findobj(bb,'Tag','Lower Adjacent Value');

            percen = 97.5;
            for i=1:num_rc
                data_sub = fc_now(:,num_rc-i+1);
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

    %         h = findobj(gca,'Tag','Box');
    %         for j=1:num_rc
    %            patch(get(h(j),'XData'),get(h(j),'YData'),col_list(2*(c-1)+ii,:),...
    %                'FaceAlpha',0.5); 
    %         end

        end
        
        if c/2<0.1
            ylim([0 2]);
        else
            ylim([0 4]);
        end
        
        ax = gca;
        hold on;
        line(linspace(ax.XLim(1),ax.XLim(2),10),ones(1,10),...
            'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
        cmb_names = {[grp_names{cmb(idx_cmb(c,1),2)} '_' grp_names{cmb(idx_cmb(c,1),1)}],...
            [grp_names{cmb(idx_cmb(c,2),2)} '_' grp_names{cmb(idx_cmb(c,2),1)}]};
        ylabel(['FC (' cmb_names{1} ' and ' cmb_names{2} ')'], 'Interpreter','none');
        set(findobj(gca,'type','line'),'linew',1);
        set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
        ax.XTickLabelRotation = 90;

        fig.PaperUnits = 'inches';
        fig.PaperSize = [num_rc*tmp+0.5 num_rc*tmp];
        fig.PaperPosition = [0 0 num_rc*tmp+0.5 num_rc*tmp];
    %     saveas( gcf, [ savedir '/v_fc_' cmb_names{1} '__' cmb_names{2} '.png' ] );
        saveas( gcf, [ savedir '/v_fc_' cmb_names{1} '__' cmb_names{2} '.pdf' ] );
        close all;

    end

end
