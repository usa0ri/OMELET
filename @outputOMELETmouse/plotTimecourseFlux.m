function plotTimecourseFlux(obj,savedir_)
   
    savedir = [savedir_ '/timecourse_flux'];
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
    
    percent_now = 97.5;
    mousetype = model_data.mousetype;
    timepoints = model_data.timepoints;
    col_list = model_data.col(idx_s,:);
    fig = figure();
    fig.PaperSize = [40 5];
    fig.PaperPosition = [0 0 40 5];
    
    hold on;
    for m=1:length(mousetype)
        idx_m = contains(grp_names,mousetype{m});
        col_now = mean(col_list(idx_m,:),1);
        for i=1:model_data.out.num_rc
            subplot(length(mousetype),model_data.out.num_rc,...
                model_data.out.num_rc*(m-1)+i);
            hold on;
            
            v_now = par.v(:,idx_m,i);
            v_median = median(v_now,1);
            v_upper = prctile(v_now,percent_now,1);
            v_lower = prctile(v_now,100-percent_now,1);
            
            ar = area(timepoints,...
                [v_lower; v_upper-v_lower]');
            set(ar(1),'FaceColor','none','EdgeColor','none');
            set(ar(2),'FaceColor',col_now,'FaceAlpha',0.1,'EdgeColor','none');
            
            line(timepoints,v_median,...
                'Color',mean(col_list(idx_m,:),1),...
                'LineWidth',2);
            scatter(timepoints,v_median,'filled','CData',col_now,'SizeData',5);
            
%             ylim([0 3.5]);
            
            xlabel('Time (h)');
            ylabel('Metabolic flux');
            title([fluxnames{i} ' (' mousetype{m} ')']);
        end
        
    end
    
    fig.PaperUnits = 'inches';
    
    set(findobj(gca,'type','axes'),'FontSize',10);
    
    print(fig,'-painters',...
        [savedir '/timecourse_flux.pdf'],...
        '-dpdf');
    
end