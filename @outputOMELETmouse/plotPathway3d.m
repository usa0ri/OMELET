function plotPathway3d(obj,bg_col)
    
    [data_layers,data_edges] = make_data_layers(obj);
    
    savedir = obj.coord3d_pathway.savedir;
    plot_data(obj,bg_col,data_layers,data_edges,savedir);
    
    % visualize the subnetwork of the target reaction
    plot_data_sub(obj,bg_col,data_layers,data_edges,savedir);
    
end


function [data_layers,data_edges] = make_data_layers(obj)
    
    % load omics data
    omics_data = obj.omics_data;
    
    % make dictionary
    met_names = obj.coord3d_pathway.met_names;
    rxn_names = obj.coord3d_pathway.rxn_names;
    [data_met, data_pro, data_rna] = make_data_omics(omics_data,...
        met_names,rxn_names);
    
    % Pyr and OAA
    idx_g = obj.model_data.idx_g;
    for g=1:obj.model_data.num_g
        data_met(ismember(met_names,'Pyruvate'),idx_g(g,1):idx_g(g,2))...
            = mean(obj.par.c_Pyr_out(:,g));
        data_met(ismember(met_names,'OAA'),idx_g(g,1):idx_g(g,2))...
            = mean(obj.par.c_OAA_out(:,g));
    end
    
    % cofactors and allosteric effectors
    met_names_eff = obj.coord3d_pathway.met_names_eff;
    [data_met_eff, ~, ~] = make_data_omics(omics_data,met_names_eff,[]);
    
    data_list{1} = [data_met;data_met_eff];
    data_list{2} = data_pro;
    data_list{3} = data_rna;
    
    % calculate fold change
    cmb = [1 3;2 4;1 2;3 4];
    fc_list = cell(1,length(data_list));
    for i=1:length(data_list)
        data_now = data_list{i};
        fc_now = nan(size(data_now,1),size(cmb,1));
        for c=1:size(cmb,1)
           fc_now(:,c) =  nanmean(data_now(:,idx_g(cmb(c,2),1):idx_g(cmb(c,2),2)),2)...
            ./nanmean(data_now(:,idx_g(cmb(c,1),1):idx_g(cmb(c,1),2)),2);
        end
        fc_list{i} = fc_now;
    end
    
    % flux
    idx_cmb = [2 5 1 6];
    fc_flux = reshape(mean(obj.fc_wtob(:,:,idx_cmb)),length(rxn_names),length(idx_cmb));
    fc_flux = fc_flux(obj.fig_info.idx_reorder,:);% FIXME
    
    % output
    data_layers{1} = fc_list{1};
    data_layers{2} = fc_flux;
    data_layers{3} = fc_list{2};
    data_layers{4} = fc_list{3};
    
    %%%%%%%%%%%%%%%%%
    num_species = size(obj.cont_flux.intgrp,2);
    cont_intgrp = nan(length(rxn_names),length(idx_cmb),num_species);
%     % scaling contribution by flux fold change
%     for s=1:num_species
%         cont_intgrp(:,:,s) = reshape(mean(obj.cont_flux.intgrp(:,s,:,idx_cmb),[2 3]),...
%             length(rxn_names),length(idx_cmb)) .* fc_flux;
%     end
    % NOT scaling
    for s=1:num_species
        cont_intgrp(:,:,s) = reshape(mean(obj.cont_flux.intgrp(:,s,:,idx_cmb),[2 3]),...
            length(rxn_names),length(idx_cmb));
    end
    %%%%%%%%%%%%%%%%%
    
    % contribution from metabolite to flux
    sto_list = fieldnames(obj.coord3d_pathway.sto);
    num_edges = length(obj.coord3d_pathway.edges);
    data_edges = cell(1,num_edges);
    % substrate/product/cofactors/allosteric effectors - flux
    for i=1:num_edges
        if i==1
%             data_edges{i} = reshape(mean(obj.cont.intgrp(:,1,:,idx_cmb),[2 3]),...
%                 length(rxn_names),length(idx_cmb));
            data_edges{i} = sum(cont_intgrp(:,:,[1 2]),3);
        elseif i==2
%             data_edges{i} = reshape(mean(obj.cont_flux.intgrp(:,1,:,idx_cmb),[2 3]),...
%                 length(rxn_names),length(idx_cmb));
            data_edges{i} = cont_intgrp(:,:,1);
        else
%             tmp = reshape(mean(obj.cont_flux.intgrp(:,i,:,idx_cmb),[2 3]),...
%                 length(rxn_names),length(idx_cmb));
            tmp = cont_intgrp(:,:,i);
            S_now = getfield(obj.coord3d_pathway.sto,sto_list{i-2});
            data_edges_ = cell(size(S_now,1),1);
            for m=1:size(S_now,1)
                idx_ = find(S_now(m,:)>0);
                data_edges_{m} = tmp(idx_,:);
            end
            data_edges{i} = data_edges_;
        end
    end
    
end

function [data_met, data_pro, data_rna] = make_data_omics(omics_data,...
    met_names,rxn_names)

%     metabolome data
    met_names_tmp = omics_data.var_omics{1};
    met_names_tmp{ismember(met_names_tmp,'G3P')} = 'Glycerol 3P';
    met_names_tmp{ismember(met_names_tmp,'Fumarate')} = 'Fum';
    met_names_tmp{ismember(met_names_tmp,'Citrate')} = 'Cit';
    met_names_tmp{ismember(met_names_tmp,'Succinate')} = 'Suc';
    met_names_tmp{ismember(met_names_tmp,'Acetyl-CoA')} = 'AcCoA';
    met_names_tmp{ismember(met_names_tmp,'Malate')} = 'Mal';
    % met_names_tmp{ismember(met_names_tmp,'NAD')} = 'NAD+';

    met_names_ = cellfun(@(x) char(extractBefore(x,'_')), met_names, 'UniformOutput',false);
    [~,idx_now] = ismember(met_names,met_names_tmp);
    [~,idx_now2] = ismember(met_names_,met_names_tmp);
    idx_now(idx_now2>0) = nonzeros(idx_now2);

    idx_m = idx_now;
    % idx_m(idx_m>0) = idx_m_(nonzeros(idx_m));

    data_m_all = omics_data.data_omics{1};
    data_m = data_m_all(:,nonzeros(idx_m));
    
    num_smpl = size(data_m,1);
    data_met = nan(length(met_names),num_smpl);
    data_met(idx_m>0,:) = data_m';
    
%     proteome and transcriptome data

    if ~isempty(rxn_names)
%         proteome
        [~,idx_p] = ismember(rxn_names,omics_data.var_omics{2});
        idx_p(ismember(rxn_names,{'Gpt','Glud1'})) = 0;
        valset_p = rxn_names;
        valset_p(idx_p==0) = cellfun(@(x) [x '_noData'],valset_p(idx_p==0),'UniformOutput',false);

        data_p_all = omics_data.data_omics{2};
        data_p = data_p_all(:,nonzeros(idx_p));
        
%         transcriptome
        [~,idx_t] = ismember(rxn_names,omics_data.var_omics{3});
        valset_t = rxn_names;
        valset_t(idx_t==0) = cellfun(@(x) [x '_noData'],valset_t(idx_t==0),'UniformOutput',false);
        
        data_t_all = omics_data.data_omics{3};
        data_t = data_t_all(:,nonzeros(idx_t));
        
        data_pro = nan(length(rxn_names),num_smpl);
        data_rna = nan(length(rxn_names),num_smpl);
        data_pro(idx_p>0,:) = data_p';
        data_rna(idx_t>0,:) = data_t';
    else
        data_pro = nan(length(rxn_names),num_smpl);
        data_rna = nan(length(rxn_names),num_smpl);
    end

end

function plot_data(obj,bg_col,data_layers,data_edges,savedir)
    
    data_xy = data_layers;
    data_edge = data_edges;
    
    out = obj.coord3d_pathway;
    node_list = out.node_list;
    edge_list = out.edge_list;
    
    cmb_names = {'WT0hob0h','WT4hob4h','WT0hWT4h','ob0hob4h'};
    
    for c=1:length(cmb_names)
        data_xy_now = cellfun(@(x) x(:,c),data_xy,'UniformOutput',false);
        data_edge_now = cell(1,length(data_edge));
        for i=1:length(data_edge)
            if i<3 % not from metabolites
                data_edge_now{i} = data_edge{i}(:,c);
            else
                data_edge_now{i} = cellfun(@(x) x(:,c),data_edge{i},'UniformOutput',false);
            end
        end
        
        fig = figure('visible','off','Color',bg_col);
        hold on;
        % edges
        lib.plot_transomics3d('edges',edge_list,fig,bg_col,...
            out.edges,data_edge_now);
        % nodes
        lib.plot_transomics3d('multiLayer',node_list,fig,bg_col,...
            out.xy_layers,out.line_layers(2,:),out.line_layers(1,:),out.frame_layers,data_xy_now);

        ax = gca;
        ax.Visible = 'off';
        ax.View = [45 45];
        savefig([savedir '/network3d_' cmb_names{c}]);
%         keyboard;
        close all;
        
    end
    
end

function plot_data_sub(obj,bg_col,data_layers,data_edges,savedir)
    
    data_xy = data_layers;
    data_edge = data_edges;
    
    out = obj.coord3d_pathway;
    node_list = out.node_list;
    edge_list = out.edge_list;
    
    rxn_names = obj.coord3d_pathway.rxn_names;
    sub_list_names = ['PyrCycling';'GNG+Pyr';'GNG';'TCA'; rxn_names];
    sub_list_rxn_ = {{'Pklr','Pck1','Pcx'};...
        {'Gpi1','Fbp1','Gpd1','Pgam1','Eno1','Ldha','Gpt'};...
        {'Gpi1','Fbp1','Gpd1','Pgam1','Eno1'};...
        {'Pgm2','Cs','Sdha','Fh1','Mdh2','Glud1'}};
    sub_list_rxn = [sub_list_rxn_; rxn_names];
    data_sub = make_data_sub(obj,data_xy,data_edge,sub_list_rxn);
    
    cmb_names = {'WT0hob0h','WT4hob4h','WT0hWT4h','ob0hob4h'};
    
    for c=1:length(cmb_names)
%         data_xy_now = cellfun(@(x) x(:,c),data_xy,'UniformOutput',false);
        for i=1:length(sub_list_rxn)
            data_edge_now = data_sub.edge_cmb{c}{i};
            data_node_now = data_sub.node_cmb{c}{i};
            fig = figure('visible','off','Color',bg_col);
            hold on;
            % edges
            lib.plot_transomics3d('edges_sub',edge_list,fig,bg_col,...
                out.edges,data_edge_now);
            % nodes
            lib.plot_transomics3d('multiLayer',node_list,fig,bg_col,...
                out.xy_layers,out.line_layers(2,:),out.line_layers(1,:),out.frame_layers,data_node_now);

            ax = gca;
            ax.Visible = 'off';
            ax.View = [45 45];
            savefig([savedir '/subnetwork3d_' sub_list_names{i} '_' cmb_names{c}]);
    %         keyboard;
            close all;       
        end
 
    end
    
end

function data_sub = make_data_sub(obj,data_xy,data_edge,sub_list_rxn)
    
    rxn_names = obj.coord3d_pathway.rxn_names;
    sto_list = fieldnames(obj.coord3d_pathway.sto);
    
    data_edge_sub_list = cell(1,length(sub_list_rxn));
    data_node_sub_list = cell(1,length(sub_list_rxn));
    for r=1:length(sub_list_rxn)
        rxn_now = sub_list_rxn{r};
        [~,idx_r] = ismember(rxn_now,rxn_names);
        data_edge_now = cell(size(data_edge));
        data_node_now = cell(size(data_xy));
        tmp_node_m = nan(size(data_xy{1}));% nodes of metabolites
        for i=1:length(data_edge_now)% for each edge
            if i==1% node data for enzyme proteins and transcripts
                for j=2:4
                    tmp_node = nan(size(data_xy{j}));
                    tmp_node(idx_r,:) = data_xy{j}(idx_r,:);
                    data_node_now{j} = tmp_node;
                end
            end
            if i<3% enzyme protein-flux, enzyme transcripts-protein
               tmp_edge = nan(size(data_edge{i}));
               tmp_edge(idx_r,:) = data_edge{i}(idx_r,:);
               data_edge_now{i} = tmp_edge;
               
            else% metabolite - flux
                tmp_edge = cell(size(data_edge{i}));
                S_now = getfield(obj.coord3d_pathway.sto,sto_list{i-2});
                % metabolites connected to the reaction
                [idx_m,~] = find(S_now(:,idx_r)~=0);
                tmp_node_m(idx_m,:) = data_xy{1}(idx_m,:);
                
                for ii=1:length(idx_m)
                    % edges connected to the metabolite
                    data_edge_m = data_edge{i}{idx_m(ii)};
                    if size(data_edge_m,1)>1
                        idx_r_ = find(S_now(idx_m(ii),:));
                        idx_ = arrayfun(@(x) x==idx_r_,idx_r,'UniformOutput',false);
                        % if the metabolite is connected to several reactions,
                        % extract data of only the target reaction
                        for iii=1:length(idx_r)
                            if any(idx_{iii})
                                tmp_edge_ = nan(size(data_edge_m));
                                tmp_edge_(idx_{iii},:) = data_edge_m(idx_{iii},:);
                                tmp_edge{idx_m(ii)} = tmp_edge_;
                            end
                        end 
                    else
                        tmp_edge{idx_m(ii)} = data_edge_m;
                    end
                end
                data_edge_now{i} = tmp_edge;
            end
        end
        data_node_now{1} = tmp_node_m;
        data_edge_sub_list{r} = data_edge_now;
        data_node_sub_list{r} = data_node_now;
    end
    
    % organize by combination
    num_cmb = size(tmp_node_m,2);
    data_edge_sub_list_cmb = cell(1,num_cmb);
    data_node_sub_list_cmb = cell(1,num_cmb);
    for c=1:num_cmb
        data_edge_now = cell(1,length(sub_list_rxn));
        data_node_now = cell(1,length(sub_list_rxn));
        for r=1:length(sub_list_rxn)
            data_edge_sub_list_now = data_edge_sub_list{r};
            tmp_edges = cell(size(data_edge_sub_list_now));
            for i=1:length(data_edge_sub_list_now)
               if i<3
                  tmp_edges{i} = data_edge_sub_list_now{i}(:,c); 
               else
                   data_edge_sub_list_now2 = data_edge_sub_list_now{i};
                   tmp_edges2 = cell(size(data_edge_sub_list_now2));
                   for j=1:length(data_edge_sub_list_now2)
                      if ~isempty(data_edge_sub_list_now2{j}) 
                          tmp_edges2{j} = data_edge_sub_list_now2{j}(:,c);
                      end
                   end
                    tmp_edges{i} = tmp_edges2;
               end
            end
            data_node_sub_list_now = data_node_sub_list{r};
            tmp_nodes = cellfun(@(x) x(:,c),data_node_sub_list_now,...
                'UniformOutput',false);
            data_edge_now{r} = tmp_edges;
            data_node_now{r} = tmp_nodes;
        end
        data_edge_sub_list_cmb{c} = data_edge_now;
        data_node_sub_list_cmb{c} = data_node_now;
    end
    
    data_sub.edge = data_edge_sub_list;
    data_sub.node = data_node_sub_list;
    data_sub.edge_cmb = data_edge_sub_list_cmb;
    data_sub.node_cmb = data_node_sub_list_cmb;
    
end

function data_node_sub_list = make_data_node_sub(obj,data_edge_sub_list)
    
    % xyz coordinates of edges
    edges = obj.coord3d_pathway.edges;
    
    data_node_sub_list = cell(size(data_edge_sub_list));
    for i=1:length(data_edge_sub_list)
        data_edge_sub_list_now = data_edge_sub_list{i};
        data_node_sub_list_now = cell(size(data_edge_sub_list_now));
        for j=1:length(data_edge_sub_list_now)
            data_edge_now = data_edge_sub_list_now{j};
            edge_now = edges{j};
            if j<3
                idx_edge = all(~isnan(data_edge_now),2);
                edge_ = edge_now(idx_edge,:);
            else
                
            end
            
        end
        data_node_sub_list{i} = data_node_sub_list_now;
    end
    
    
    
end
