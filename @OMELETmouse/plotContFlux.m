function plotContFlux(obj,savedir)
   
    mkdir(savedir);
    
    % scatter plot of contribution - flux fold changes

    % flux fold changes
    grp_names = obj.model_data.grp_names;
    fluxnames = obj.model_data.out.rxn_names_include;
    idx_cmb = [2 5 1 6];
    cmb = obj.cont.cmb([2 5 1 6],:);
    cmb_names = {'WT0h_ob0h','WT4h_ob4h','WT0h_WT4h','ob0h_0b4h'};
    iter = size(obj.par.v,1);
    v_fc = nan(iter,length(cmb_names),length(fluxnames));
    for c=1:length(cmb_names)
        v_fc(:,c,:) = obj.par.v(:,cmb(c,2),:) ./...
            obj.par.v(:,cmb(c,1),:);
    end
    mean_v_fc = reshape(mean(v_fc,1),length(cmb_names),length(fluxnames));
    std_v_fc = reshape(std(v_fc,[],1),length(cmb_names),length(fluxnames));
    
    % contributions
    % transcripts, proteins, metabolites
    reac_names = {'transcript','protein','metabolite'};
    idx_reactant = {1, 1:2, 3:10};
    cont_flux_mat = nan(length(fluxnames),length(reac_names),iter,length(cmb_names));
    for r=1:length(reac_names)
        cont_flux_mat(:,r,:,:) = nansum(obj.cont_flux.intgrp(:,idx_reactant{r},:,idx_cmb),2);
    end
    mean_cont_flux_mat = reshape(mean(cont_flux_mat,3),...
        length(fluxnames),length(reac_names),length(cmb_names));
    std_cont_flux_mat = reshape(std(cont_flux_mat,[],3),...
        length(fluxnames),length(reac_names),length(cmb_names));
    col_reactant = [lib.my_colors('matlab_lblue',1,false);...
        lib.my_colors('matlab_blue',1,false);...
        1 0 1];
    
    % plot all MCMC samples
%     num_sub = ceil(sqrt(length(fluxnames)));
%     for c=1:length(cmb_names)
%         for r=1:length(reac_names)
%             fig = figure();
%             hold on;
%             fig.PaperSize = [10 10];
%             ax = gca;
%             lib.set_minmargin(fig,ax);
%             for i=1:length(fluxnames)
%                 subplot(num_sub,num_sub,i);
%                 v_fc_now = v_fc(:,c,i);
%                 cont_flux_mat_now = cont_flux_mat(i,r,:,c);
%                 scatter(v_fc_now(:),cont_flux_mat_now(:),...
%                     'SizeData',2,...
%                     'MarkerFaceColor',col_reactant(r,:),...
%                     'MarkerEdgeColor',col_reactant(r,:));
%                 hold on;
%                 scatter(mean_v_fc(c,i),mean_cont_flux_mat(i,r,c),...
%                     'SizeData',2,...
%                     'MarkerFaceColor',[0 0 0],...
%                     'MarkerEdgeColor',[0 0 0]);
%                 
%                 title(fluxnames{i});
%                 xlabel('flux fold change');
%                 ylabel('contribution');
%                 xlim([0 5]);
%                 ylim([0 1]);
%                 set(findobj(gca,'type','axes'),'FontSize',5,'FontName','SansSerif');
%             end
%             
%             sgtitle(['Contribution of ' reac_names{r} ' (' cmb_names{c}  ')'],'Interpreter','none');
%             print(fig,'-painters',...
%                 [savedir '/scatter_mcmc_' reac_names{r} '_' cmb_names{c} '.pdf'],...
%                 '-dpdf','-bestfit');
%             close all;
%         end 
%     end
    
    % plot mean
    for c=1:length(cmb_names)
        fig = figure('visible','off');
        fig.PaperSize = [20 20];
        ax = gca;
        lib.set_minmargin(fig,ax);
        hold on;
        for r=1:length(reac_names)
            subplot(2,2,r);
            v_now = mean_v_fc(c,:);
            cont_now = mean_cont_flux_mat(:,r,c);
            hold on;
            line(repmat(1.5,1,2),[0 1],...
                'LineWidth',1,'Color','k');
            eb(1) = errorbar(v_now,cont_now,...
                std_v_fc(c,:), 'horizontal', 'LineStyle', 'none',...
                'Color',[0.6 0.6 0.6],'LineWidth',1);
            eb(2) = errorbar(v_now,cont_now,...
                std_cont_flux_mat(:,r,c), 'vertical', 'LineStyle', 'none',...
                'Color',[0.6 0.6 0.6],'LineWidth',1);
            scatter(v_now,cont_now,...
                'SizeData',10,...
                'MarkerFaceColor',col_reactant(r,:),...
                'MarkerEdgeColor',col_reactant(r,:));
            
            text(v_now,cont_now,fluxnames,'FontSize',8,'Color','k');
            title(['Contribution of ' reac_names{r}]);
            xlabel('Flux fold change');
            ylabel('Contribution');
            xlim([0 3.5]);
            xticks(0:0.5:3.5);
            ylim([0 1]);
            set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');
        end
        sgtitle(cmb_names{c},'Interpreter','none');
        print(fig,'-painters',...
            [savedir '/scatter_' cmb_names{c} '.pdf'],...
            '-dpdf','-bestfit');
        close all;
    end
    
end