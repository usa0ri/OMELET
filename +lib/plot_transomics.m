function plot_transomics(type,name,fig,varargin)
    
    switch type
        case 'singleLayer'
%             xy,line,line_all,frame,data_xy
            plot_single(name,fig,varargin);
            
        case 'multiLayer'
%             name = {'met','flux','pro','rna'};
            num_layer = length(varargin{1});
            for i=1:num_layer
                args_now = cellfun(@(x) x{i}, varargin, 'UniformOutput',false);
                plot_single(name{i},fig,args_now);
            end
            
        case 'edges'
            num_edges = length(varargin{1});
            for i=1:num_edges
                args_now = cellfun(@(x) x{i}, varargin, 'UniformOutput',false);
                plot_edge(name{i},fig,args_now);
            end
    end
    
    
end

function plot_single(node_name,fig,args)
    
    n = 2;
    xy = args{1};
    line_rxn = args{2};
    line_rxn_all = args{3};
    frame = args{4};
    data_xy = args{5};
    
    figure(fig);
    hold on;
    
    switch node_name
        case 'met'
            xy_col = [1 0 0];
        case 'flux'
            xy_col = [];
        case 'pro'
            xy_col = [0 0 1];
        case 'rna'
            xy_col = [0 0 1];
    end
    
    col_line_all = [0.3 0.3 0.3];
    for r=1:size(line_rxn_all,1)
        line(linspace(line_rxn_all(r,1),line_rxn_all(r,3),n),...
            linspace(line_rxn_all(r,2),line_rxn_all(r,4),n),...
            'Color',col_line_all,'LineWidth',1);
    end
    
    % if flux
    col_flux = [0.5 0.5 0.5];
    if isempty(xy_col) && ~all(data_xy==0)
        for r=1:size(line_rxn,1)
            line(linspace(line_rxn(r,1),line_rxn(r,3),n),...
                linspace(line_rxn(r,2),line_rxn(r,4),n),...
                'Color',col_flux,'LineWidth',data_xy(r)*2.5);
        end
    else
        for r=1:size(line_rxn,1)
            line(linspace(line_rxn(r,1),line_rxn(r,3),n),...
                linspace(line_rxn(r,2),line_rxn(r,4),n),...
                'Color',col_flux,'LineWidth',1);
        end
    end
        
    if ~isempty(xy) && ~isempty(xy_col)
        if isempty(data_xy)
            scatter(xy(:,1),xy(:,2),'filled','MarkerFaceColor',xy_col,...
            'SizeData',10);
        else
            for i=1:size(xy,1)
                scatter(xy(i,1),xy(i,2),'filled','MarkerFaceColor',xy_col,...
                    'SizeData',data_xy(i)*20);
            end
            
        end
        
    end
    for f=1:4
        ff = f+1;
        if f==4
           ff = 1; 
        end
        line(linspace(frame(f,1),frame(ff,1),n),...
            linspace(frame(f,2),frame(ff,2),n),...
            'Color','k','LineWidth',1);
    end
    
end

function plot_edge(edge_name,fig,args)
    
    n=2;
    edges = args{1};
    data_edges = args{2};
    
    figure(fig);
    hold on;

    addpath ../../myMatlab/;
    switch edge_name
        case 's-flux'
            edge_col = m.my_colors('matlab_orange',1,false);
        case 'p-flux'
            edge_col = m.my_colors('matlab_yellow',1,false);
        case 'pro-flux'
            edge_col = m.my_colors('matlab_blue',1,false);
        case 'pro-rna'
            edge_col = m.my_colors('matlab_lblue',1,false);
        case 'co-flux'
            edge_col = m.my_colors('matlab_red',1,false);
        case 'allo-flux'
            edge_col = m.my_colors('matlab_purple',1,false);
    end
    
    if iscell(edges)
        for i=1:length(edges)
            edge_now = edges{i};
            if isempty(data_edges)
                data_now = [];
            else
                data_now = data_edges(i);
            end
            if ~isempty(edge_now)
                if isempty(data_now) || all(data_now==0)
                    for e=1:size(edge_now,1)
                        line(linspace(edge_now(e,1),edge_now(e,3),n),...
                            linspace(edge_now(e,2),edge_now(e,4),n),...
                            'Color',edge_col,...
                            'LineWidth',1);
                    end
                else
                    idx_nonzero = find(data_now~=0);
                    for j=1:length(idx_nonzero)
                        e = idx_nonzero(j);
                        line(linspace(edge_now(e,1),edge_now(e,3),n),...
                            linspace(edge_now(e,2),edge_now(e,4),n),...
                            'Color',edge_col,...
                            'LineWidth',data_now(e)*8);
                    end
                end
            end
        end
    else
        edge_now = edges;
        data_now = data_edges;
        if isempty(data_edges)
            for e=1:size(edge_now,1)
                line(linspace(edge_now(e,1),edge_now(e,3),n),...
                    linspace(edge_now(e,2),edge_now(e,4),n),...
                    'Color',edge_col,...
                    'LineWidth',1);
            end
        else
            for e=1:size(edge_now,1)
                line(linspace(edge_now(e,1),edge_now(e,3),n),...
                    linspace(edge_now(e,2),edge_now(e,4),n),...
                    'Color',edge_col,...
                    'LineWidth',data_now(e)*8);
            end
        end
        
    end
    
end