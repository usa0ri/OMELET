function plotContViolin(obj,type,savedir)

model_data = obj.model_data;

mkdir(savedir);

grp_names = model_data.grp_names;
num_g = length(grp_names);
cmb = obj.cont.cmb;

opts = load_s_names(obj,type);

% fname_str = [savedir_now '/total'];
% plot_cont(opts.cont.total,fname_str,model_data,opts);
% 
% % group
% for g=1:num_g
%     fname_str = [savedir_now '/' grp_names{g} ];
%     plot_cont(opts.cont.grp(:,:,:,g),fname_str,model_data,opts);
% end

% inter-group
for c=1:size(cmb,1)
    fname_str = [savedir '/' grp_names{cmb(c,1)} '_' grp_names{cmb(c,2)} ];
    plot_cont(opts.cont.intgrp(:,:,:,c),fname_str,model_data,opts);
end


end


function plot_cont(cont_r_list,fname_str,model_data,opts)

    fluxnames = model_data.X.rxn.rxn_names_include;
    num_rc = model_data.X.num.num_rc;
    num_r = length(opts.s_names);
    sz = size(cont_r_list);
    cont = nan(sz(1),num_r,sz(3));
    for i=1:num_r
       cont(:,i,:) = nansum(cont_r_list(:,opts.idx_s{i},:),2); 
    end
    iter = size(cont,3);
    
    %%%%%%%%%%%%
    % reorder fluxnames
    fluxnames = fluxnames(opts.idx_reorder);
    cont = cont(opts.idx_reorder,:);
    %%%%%%%%%%%%

    fig = figure('visible','off');
    for i=1:num_rc
        sub = subplot(num_rc,1,i);
        cont_now = reshape(cont(i,:,:),num_r,iter);
        [h,L,MX,MED,bw]=lib.my_violin(cont_now',...
            'xlabel',opts.s_names,...
            'facecolor',opts.col,...
            'edgecolor',opts.col,...
            'medc','k',...
            'facealpha',1);
        ax = gca;
%         for j=1:num_r
%             x_now = ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))/(num_r+2)*(i+1);
%             text(x_now,MED(j),num2str(round(MED(j),3)),...
%                 'Color','r','FontSize',5);
%         end
        
        if i<num_rc
           ax.XTickLabel = ''; 
        end
        ax.YLim = [-0.2 1.2];
        ax.XLim = ax.XLim + 0.025;
%         xpos = linspace(ax.XLim(1),ax.XLim(2),10);
%         line(xpos,ones(1,10),'Color','k');
%         line(xpos,zeros(1,10),'Color','k');
        ax.XTickLabelRotation = 30;
        ax.YLabel.String = fluxnames{i};
    end
    fig.PaperUnits = 'inches';
    fig.PaperSize = [4 num_rc];
    fig.PaperPosition = [0 0 4 num_rc];
%     saveas( gcf, [ fname_str '.png' ] );
    saveas( gcf, [ fname_str '.pdf' ] );
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
       s_names = {'Enzyme(transcript)','Enzyme(unaccounted)',...
           'Substrate','Product','Cofactor','Allosteric effectors','Unaccounted'}; 
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