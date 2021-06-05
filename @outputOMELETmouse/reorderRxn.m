function reorderRxn(obj)
   
    if isempty(obj.fig_info)
        fluxname_new_ = {'Pgm2','Gpi1','Fbp1','Gpd1','Pgam1','Eno1',...
            'Ldha','Gpt',...
            'Pklr','Pcx','Pck1',...
            'Cs','Sdha','Fh1','Mdh2','Glud1'}';
        fluxname_new = {'PGM','GPI','FBPase','GPD','PGAM','ENO',...
            'LDH','GPT',...
            'PK','PC','PEPCK',...
            'CS','SDH','FH','MDH','GLUD'}';
        fig_info.fluxnames = fluxname_new;
    else
        fluxname_new = obj.fig_info.fluxnames;
        fig_info.fluxnames = fluxname_new;
    end
    fluxname_old = obj.model_data.out.rxn_names_include;
    
    [~,idx_] = ismember(fluxname_new_,fluxname_old);
%     fluxname_old(idx_)

    fig_info.idx_reorder = idx_;
    obj.fig_info = fig_info;
    
end