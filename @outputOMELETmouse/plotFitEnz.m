function plotFitEnz(obj,savedir)

    model_data = obj.model_data;

    savedir_now = savedir;
    mkdir(savedir_now);

    rxn_names = model_data.X.rxn.rxn_names_include;
    enz = model_data.out.enz;
    enz(ismember(rxn_names,obj.model_data.out.enz_names_est),:) = nan;

    rna = model_data.out.rna;
    rna_eff = model_data.out.rna_eff;

    num_g = model_data.num_g;
    idx_g = model_data.idx_g;
    num_smpl = model_data.num_smpl;
    num_smpl_g = model_data.num_smpl_g;
    iter = size(obj.par.a,1);
%     rna_eff_pred = nan(iter,1,num_smpl);
%     for i=1:iter
%         rna_eff_pred(i,1,:) = randn(1,num_smpl) * sqrt(obj.par.sigma_n2(i)) + obj.par.y_eff(i,:);
%     end
%     

    fnames = {'protein','transcript'};
    % is_mean = {'mean','raw'};
    is_mean = {'raw'};
    for i=1:2
       if i==1
            x_pred = obj.par.enz_pred;
            x = enz;
            % reorder fluxnames
            idx_reorder = obj.fig_info.idx_reorder;
       elseif i==2
    %         x_pred = obj.par.enz_pred2;
    %         x = enz;
    %         idx_reorder = obj.fig_info.idx_reorder;
%             x_pred = cat(2,obj.par.rna_pred,rna_eff_pred);
%             x = [rna; rna_eff];
            x_pred = obj.par.rna_pred;
            x = rna;
            idx_reorder = obj.fig_info.idx_reorder;
       end
       plot_box(obj,x_pred,model_data,x,fnames{i},is_mean{1},idx_reorder,...
               savedir_now);

    end


end

function plot_box(obj,x_pred,model_data,x,fname,is_mean,idx_,savedir)

    rxn_names = model_data.X.rxn.rxn_names_include;
    num_rc = model_data.X.num.num_include;
    iter = size(x_pred,1);
    num_smpl = model_data.out.num_smpl;
    idx_g = model_data.out.idx_g;
    num_g = model_data.out.num_g;

    %%%%%%%%%%%
    % reorder fluxnames
    switch fname
        case 'protein'
            rxn_names = rxn_names(idx_);
        case 'transcript'
            rxn_names = rxn_names(idx_);
%             rxn_names_tmp = [rxn_names; 'Sdhb'];
%             rxn_names = rxn_names_tmp(idx_);
            num_rc = num_rc;
    end
    x_pred = x_pred(:,idx_,:);
    x = x(idx_,:);
    %%%%%%%%%%%

    switch is_mean
        case 'mean'
            x_vec = nan(iter,num_rc,num_g);
            for i=1:num_rc
                x_now = x_pred(:,i,:);
                for g=1:num_g
                    x_vec(:,i,g) = nanmean(x_now(:,1,idx_g(g,1):idx_g(g,2)),3);
                end
            end
        case 'raw'
            x_vec = nan(iter*12,num_rc,num_g);
            for i=1:num_rc
                x_now = x_pred(:,i,:);
                for g=1:num_g
                    tmp = x_now(:,idx_g(g,1):idx_g(g,2));
                    x_vec(1:iter*size(tmp,2),i,g) = tmp(:);
                end
            end
    end

    col_list = obj.model_data.col;
    grp_names = obj.model_data.grp_names;

    for g=1:num_g
        x_now = x_vec(:,:,g);
        x_data_now = x(:,idx_g(g,1):idx_g(g,2));

        fig = figure('visible','off');
        boxplot(x_now,...
            'labels',rxn_names,...
            'Colors',[0.2 0.2 0.2],...
            'symbol','',...%'OutlierSize',3,...
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
            data_sub = x_now(:,num_rc-i+1);
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

        hold on;
        for r=1:num_rc
            scatter(repmat(r,idx_g(g,2)-idx_g(g,1)+1,1),x_data_now(r,:),...
                20,'filled','MarkerFaceColor',col_list(g,:));
        end
    %     ylim([0 2.5]);
        title(grp_names{g});
        set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');

        fig.PaperUnits = 'inches';

        switch fname
            case 'protein'
                fig.PaperSize = [10 4];
                fig.PaperPosition = [0 0 10 4];
            case 'transcript'
                fig.PaperSize = [10.5 4];
                fig.PaperPosition = [0 0 10.5 4];
        end
        ylim([0 2.5]);

    %     saveas( gcf, [ savedir '/fit_' tmp '_' is_mean '_' fname '_' grp_names{g} '.png' ] );
        saveas( gcf, [ savedir '/fit_' is_mean '_' fname '_' grp_names{g} '.pdf' ] );
        close all;
    end


end