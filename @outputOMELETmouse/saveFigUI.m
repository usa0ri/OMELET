function saveFigUI(obj,savedir)
    
    loaddir = obj.coord3d_pathway.savedir;
    
    angle_fixed = [103 15];
    save_fig(loaddir,[],angle_fixed,savedir);
    
%     node_size = [0.5 1.5 1 2];
%     for n=1:length(node_size)
%         save_fig(loaddir,['all_empty_' num2str(node_size(n))],angle_fixed,savedir);
%     end
%     edge_size = [0.25 0.5 0.75 1];
%     for n=1:length(edge_size)
%         save_fig(loaddir,['edges_empty_' num2str(edge_size(n))],angle_fixed,savedir);
%     end
    
    save_fig(loaddir,'Pklr',angle_fixed,savedir);
%     save_fig(loaddir,'Gpi1',angle_fixed,savedir);
%     save_fig(loaddir,'Pgm2',angle_fixed,savedir);
    
    save_fig(loaddir,'PyrCycling',angle_fixed,savedir);
%     save_fig(loaddir,'GNG+Pyr',angle_fixed,savedir);
    save_fig(loaddir,'GNG',angle_fixed,savedir);
%     save_fig(loaddir,'TCA',[100 15],savedir);
    
end

function save_fig(loaddir,rxn_now,angle_now,savedir)
    
    [fig1,ax1,fname1] = open_fig(loaddir,rxn_now,'WT0hWT4h');
    pause(5);
    [fig2,ax2,fname2] = open_fig(loaddir,rxn_now,'ob0hob4h');
    ax1.View = angle_now;
    ax2.View = angle_now;
%     keyboard;
    print(fig1,'-painters',[savedir '/' fname1 '_ui'],'-dpdf','-bestfit');
    print(fig2,'-painters',[savedir '/' fname2 '_ui'],'-dpdf','-bestfit');
    close all;
    
    [fig1,ax1,fname1] = open_fig(loaddir,rxn_now,'WT0hob0h');
    pause(5);
    [fig2,ax2,fname2] = open_fig(loaddir,rxn_now,'WT4hob4h');
    ax1.View = angle_now;
    ax2.View = angle_now;
    print(fig1,'-painters',[savedir '/' fname1 '_ui'],'-dpdf','-bestfit');
    print(fig2,'-painters',[savedir '/' fname2 '_ui'],'-dpdf','-bestfit');
    close all;
    
end

function [fig,ax,fname] = open_fig(loaddir,rxn_now,cmb_now)
    
    if isempty(rxn_now)
        fname = [ 'network3d_' cmb_now];
    elseif contains(rxn_now,{'edges','all'})
        fname = rxn_now;
    else
        fname = ['subnetwork3d_' rxn_now '_' cmb_now];
    end
    fig = openfig([loaddir '/' fname '.fig'],'invisible');
    ax = gca;
    lib.set_minmargin(fig,ax);
    set(gcf,'InvertHardCopy','off');
    
end