function makePathway3d(obj,inter,bg_col,savedir)
    
    savedir = [savedir '/pathway3d_inter' num2str(inter) '_' bg_col];
    if exist(savedir,'dir')==7
        if ischar(savedir(end))
            savedir =  [savedir '1'];
        elseif isnumeric(savedir(end))
            savedir = [savedir num2str(str2double(savedir(end))+1)];
        end
    end
    mkdir(savedir);
    
    % xyz coordinates of metabolites
    out = make_coord3d_met(obj);
    
    % lines representing reactions 
    % (only reactions whose posterior distribution is calculated)
    out = make_line3d_rxn(obj,out);

    % lines representing all the reactions 
    out = make_line3d_rxn_all(obj,out);
    
    out = plot_pathway_all3d_empty(out,inter,bg_col,savedir);
    
    % edges between layers
    out = make_edges3d(obj,out);
    out = plot_edges3d_empty(out,bg_col,savedir);
    
    out.savedir = savedir;
    
    obj.coord3d_pathway = out;
    
end

function out = make_coord3d_met(obj)
    
    met_names_ = obj.model_data.X.met.met_names_all;
    met_names = [met_names_(1:20); met_names_(23:26); met_names_(21)];
    met_names_eff_ = obj.model_data.out.met_names_eff;
    met_names_eff = [met_names_eff_(3:4);...
        met_names_eff_(1:2);...
        met_names_eff_(5:7);...
        met_names_eff_(13:14)];
    
    met_names_tmp = [met_names;met_names_eff];
    met_names_tmp{ismember(met_names_tmp,'G3P')} = 'GAP';
    met_names_tmp{ismember(met_names_tmp,'Glycerol 3P')} = 'G3P';
    met_names_tmp{ismember(met_names_tmp,'Fum')} = 'Fumarate';
    met_names_tmp{ismember(met_names_tmp,'Cit')} = 'Citrate';
    met_names_tmp{ismember(met_names_tmp,'Suc')} = 'Succinate';
    met_names_tmp{ismember(met_names_tmp,'AcCoA')} = 'Acetyl-CoA';
    met_names_tmp{ismember(met_names_tmp,'Mal')} = 'Malate';
    met_names_tmp{ismember(met_names_tmp,'NAD+')} = 'NAD';
    
    
    % import coordinates of svg file
    tbl_ = readtable([obj.data_path '/pathway_matlab1.csv']);
    [~,idx_] = ismember(met_names_tmp,tbl_.Var1);
    assert(all(idx_));
    
    % scaling
    max_y = max(tbl_.Var3);
    xy_met_ = [tbl_.Var2(idx_) max_y-tbl_.Var3(idx_)];
    xy_met = ( xy_met_ - [min(xy_met_(:,1)) min(xy_met_(:,2))] )./10;
    

    
    frame = [min(xy_met(:,1))-0.2 min(xy_met(:,2))-0.2;...
        max(xy_met(:,1))+0.2 min(xy_met(:,2))-0.2;...
        max(xy_met(:,1))+0.2 max(xy_met(:,2))+0.2;...
        min(xy_met(:,1))-0.2 max(xy_met(:,2))+0.2];
    
    out.met_names = met_names;
    out.met_names_eff = met_names_eff;
    out.idx_met = 1:length(met_names);
    out.idx_met_eff = length(met_names)+1:size(xy_met,1);
    out.xy_met = [xy_met zeros(size(xy_met,1),1)];
    out.frame = [frame zeros(size(frame,1),1)];
    
end

function figID = plot_grid(figID,axID)
    
    figure(figID,'visible','off');
    
    grid_x = linspace(axID.XLim(1),axID.XLim(2),40);
    grid_y = linspace(axID.YLim(1),axID.YLim(2),40);

    for i=1:length(grid_x)
        line(repmat(grid_x(i),1,2),grid_y([1 end]));
    end
    for i=1:length(grid_y)
       line(grid_x([1 end]),repmat(grid_y(i),1,2)) 
    end
    
end

function out = make_line3d_rxn(obj,out)
    
    met_names = out.met_names;
    xy_met = out.xy_met;
    
    rxn_names = obj.model_data.out.rxn_names_include;
    
    [~,idx_m] = ismember(met_names,obj.model_data.X.met.met_names_all);
    [~,idx_r] = ismember(rxn_names,obj.model_data.X.rxn.rxn_names_all);
    S = obj.model_data.X.sto.S_org(idx_m,idx_r);
    % reactions whose fluxes are estimated
    xy_rxn = nan(length(rxn_names),3);
    line_rxn = nan(length(rxn_names),6);
    for r=1:length(rxn_names)
        idx_ = find(S(:,r)~=0);
        if r==12
            [~,idx_] = ismember({'OAA','Cit'},met_names);% Cs
        elseif r==16
            [~,idx_] = ismember({'AKG','Glu'},met_names);% Glud1
        elseif length(idx_)~=2
            error('Substrate and products cannot be determined');
        end
        xy_rxn(r,:) = mean(xy_met(idx_,:),1);
        line_rxn(r,:) = [xy_met(idx_(1),:) xy_met(idx_(2),:)];
    end
    
    out.rxn_names = rxn_names;
    out.xy_rxn = xy_rxn;
    out.line_rxn = line_rxn;
    out.S = S;

end

function out = make_line3d_rxn_all(obj,out)
    
    met_names = out.met_names;
    xy_met = out.xy_met;
    
    rxn_names_all = obj.model_data.X.rxn.rxn_names_all;
    idx_exclude = ismember(rxn_names_all,...
        {'TCAflux','Got1','Got2'});
    rxn_names_all_ = rxn_names_all(~idx_exclude);
    met_names_all = obj.model_data.X.met.met_names_all;
    S_all = obj.model_data.X.sto.S_org(:,~idx_exclude);
    line_rxn_all = nan(length(rxn_names_all_)+2,6);
    rr = 1;
    for r=1:length(rxn_names_all_)
        met_name_now = met_names_all(S_all(:,r)~=0);
        if r==17
            met_name_now = {'Cit','AcCoA','OAA'};
        elseif r==23 % Glud1
            met_name_now = {'AKG','Glu'};
        end
        [~,idx_now] = ismember(met_name_now,met_names);
        if r==6 || r==17% Aldob or Cs
            rxn_names_all = [rxn_names_all(1:r) ;...
                [rxn_names_all{r} '*'];...
                rxn_names_all((r+1):end)];
            line_rxn_all(rr,:) = [xy_met(idx_now(1),:) xy_met(idx_now(2),:)];
            rr = rr+1;
            line_rxn_all(rr,:) = [xy_met(idx_now(1),:) xy_met(idx_now(3),:)];
            rr = rr+1;
        else
            line_rxn_all(rr,:) = [xy_met(idx_now(1),:) xy_met(idx_now(2),:)];
            rr = rr+1;
        end
    end
    assert(rr==size(line_rxn_all,1)+1);
    
    out.line_rxn_all = line_rxn_all;
    out.rxn_names_all = rxn_names_all;
    out.S_all = S_all;
    
end

function out = plot_pathway_all3d_empty(out,inter,bg_col,savedir)
    
    xy_met = out.xy_met;
    xy_rxn = out.xy_rxn;
    line_rxn = out.line_rxn;
    line_rxn_all = out.line_rxn_all;
    frame = out.frame;

    % make frame lines
    frame_layer_ = [frame;frame(1:2,:)];
    
    % coordinates in each layer
    xy_layers = cell(1,4);
    frame_layers = cell(1,4);
    % lines (included and all)
    line_layers = cell(2,4);
    for i=1:4
        if i==1
            xy_layers{i} = xy_met +...
                repmat([0 0 1],size(xy_met,1),1).*inter.*(i-1);
        else
            xy_layers{i} = xy_rxn +...
                repmat([0 0 1],size(xy_rxn,1),1).*inter.*(i-1);
        end
        frame_layers{i} = frame_layer_ + ...
            repmat([0 0 1],size(frame_layer_,1),1).*inter.*(i-1); 
        line_layers{1,i} = line_rxn_all + ...
            repmat([0 0 1 0 0 1],size(line_rxn_all,1),1).*inter.*(i-1);
        line_layers{2,i} = line_rxn + ...
            repmat([0 0 1 0 0 1],size(line_rxn,1),1).*inter.*(i-1);
    end
    
    out.xy_layers = xy_layers;
    out.line_layers = line_layers;
    out.frame_layers = frame_layers;
    
    % plot
    node_size = [0.5, 1, 1.5, 2];
    for n=1:length(node_size)
        data_xy = cell(1,length(xy_layers));
        for i=1:length(xy_layers)
            data_xy{i} = ones(size(xy_layers{i},1),1).*node_size(n);
        end
        name = {'met','flux','pro','rna'};
        fig = figure('Color',bg_col,'visible','off');
        hold on;
        lib.plot_transomics3d('multiLayer',name,fig,bg_col,...
            xy_layers,line_layers(2,:),line_layers(1,:),frame_layers,data_xy);
        ax = gca;
        ax.Color = bg_col;
        ax.Visible = 'off';
        ax.View = [45 45];
        savefig([savedir '/all_empty_' num2str(node_size(n)) '.fig' ]);
        close all;
    end

end

function out = make_edges3d(obj,out)
    
    xy_layers = out.xy_layers;
    edges = cell(1,2);% FIXME
    
    % flux - enzyme protein
    % enzyme protein - transcripts
    for i=1:2
        edges{i} = [xy_layers{i+1} xy_layers{i+2}];
    end
    
    met_names_all = [out.met_names;out.met_names_eff];
    xy_rxn = xy_layers{2};
    % metabolite - flux
    % cofactors
    % stoichiometry matrix of cofactor
    m = [out.met_names;out.met_names_eff];
    r = out.rxn_names;
    S_eff = sparse(length(m),length(r));
    S_cofactor_sub = S_eff;
    S_cofactor_pro = S_eff;
    S_cofactor_sub(ismember(m,'NAD+'),ismember(r,'Gpd1')) = 1;
    S_cofactor_pro(ismember(m,'NADH'),ismember(r,'Gpd1')) = 1;
    S_cofactor_sub(ismember(m,'NAD+'),ismember(r,'Ldha')) = 1;
    S_cofactor_pro(ismember(m,'NADH'),ismember(r,'Ldha')) = 1;
    S_cofactor_sub(ismember(m,'NAD+'),ismember(r,'Mdh2')) = 1;
    S_cofactor_pro(ismember(m,'NADH'),ismember(r,'Mdh2')) = 1;
    S_cofactor_sub(ismember(m,'NAD+'),ismember(r,'Glud1')) = 1;
    S_cofactor_pro(ismember(m,'NADH'),ismember(r,'Glud1')) = 1;
    S_cofactor_sub(ismember(m,'ADP'),ismember(r,'Pklr')) = 1;
    S_cofactor_sub(ismember(m,'ATP'),ismember(r,'Pcx')) = 1;
    S_cofactor_sub(ismember(m,'GTP'),ismember(r,'Pck1')) = 1;
    S_cofactor_sub(ismember(m,'FAD'),ismember(r,'Sdha')) = 1;
    
    % allosteric effectors
    % stoichiometry matrix of allosteric effectors
    S_allo_i = S_eff;
    S_allo_a = S_eff;
    S_allo_a(ismember(m,'Cit'),ismember(r,'Fbp1')) = 1;
    S_allo_i(ismember(m,'AMP'),ismember(r,'Fbp1')) = 1;
    S_allo_a(ismember(m,'F1,6P'),ismember(r,'Pklr')) = 1;
    S_allo_i(ismember(m,'ATP'),ismember(r,'Pklr')) = 1;
    S_allo_a(ismember(m,'AcCoA'),ismember(r,'Pcx')) = 1;
    S_allo_a(ismember(m,'ADP'),ismember(r,'Glud1')) = 1;
    S_allo_i(ismember(m,'ATP'),ismember(r,'Glud1')) = 1;
    S_allo_i2 = S_eff;
    S_allo_i2(ismember(m,'Ala'),ismember(r,'Pklr')) = 1;
    S_allo_i2(ismember(m,'GTP'),ismember(r,'Glud1')) = 1;
    S_allo_i3 = S_eff;
    S_allo_i3(ismember(m,'Phe'),ismember(r,'Pklr')) = 1;
    S_allo_i3(ismember(m,'Leu'),ismember(r,'Glud1')) = 1;

    sto.S_sub = [out.S<0; zeros(length(out.met_names_eff),length(out.rxn_names))];
    sto.S_pro = [out.S>0; zeros(length(out.met_names_eff),length(out.rxn_names))];
    % irreversible
    rxn_irrev = obj.model_data.out.rxn_names(obj.model_data.out.is_irrev);
    [~,idx_irrev] = ismember([rxn_irrev;'Glud1'],out.rxn_names);
    sto.S_pro(:,idx_irrev) = 0;
    sto.S_cofactor_sub = S_cofactor_sub;
    sto.S_cofactor_pro = S_cofactor_pro;
    sto.S_allo_a = S_allo_a;
    sto.S_allo_i = S_allo_i;
    sto.S_allo_i2 = S_allo_i2;
    sto.S_allo_i3 = S_allo_i3;
    
    clear m r;
    
    sto_list = fieldnames(sto);
    edges_eff = cell(1,length(sto_list));
    for s=1:length(sto_list)
        edges_eff_now = cell(length(met_names_all),1);
        S_ = getfield(sto,sto_list{s});
        for m=1:length(met_names_all)
            idx_tmp = find(S_(m,:)~=0);
            if ~isempty(idx_tmp)
                edges_eff_ = [];
                for i=1:length(idx_tmp)
                    edges_eff_ = [edges_eff_;out.xy_met(m,:) xy_rxn(idx_tmp(i),:)];
                end
                edges_eff_now{m} = edges_eff_;
            end
        end
        edges_eff{s} = edges_eff_now;
    end
    
    edges = [edges,edges_eff];
    
    out.edges = edges;
    out.sto = sto;
    
end

function out = plot_edges3d_empty(out,bg_col,savedir)
    
    node_list = {'met','flux','pro','rna'};
    edge_list = {'pro-flux','pro-rna','s-flux','p-flux',...
        'co-flux','co-flux','allo-flux','allo-flux','allo-flux','allo-flux'};
    data_xy = cell(1,length(node_list));
    
    edge_size = [0.25, 0.5, 0.75, 1];
    for n=1:length(edge_size)
        data_edge = cell(1,length(edge_list));
        for i=1:length(edge_list)
            data_edge{i} = ones(size(out.edges{i},1),1).*edge_size(n);
        end
        
        fig = figure('Color',bg_col,'visible','off');
        hold on;
        % edges
        lib.plot_transomics3d('edges',edge_list,fig,bg_col,...
            out.edges,data_edge);
        % nodes
        lib.plot_transomics3d('multiLayer',node_list,fig,bg_col,...
            out.xy_layers,out.line_layers(2,:),out.line_layers(1,:),out.frame_layers,data_xy);
        ax = gca;
        ax.Visible = 'off';
        ax.View = [45 45];
        savefig([savedir '/edges_empty_' num2str(edge_size(n)) '.fig' ]);

        close all;
    end
    
    out.node_list = node_list;
    out.edge_list = edge_list;
    
end