function plot_transomics3d(type,name,fig,bg_col,varargin)
    
    switch type
        case 'edge'
            % xyz coordinates of edges: varargin{1}
            % data of edges: varargin{2}
            for i=1:length(varargin)
                assert(length(name)==length(varargin{i}));
            end
        case 'multiLayer'
    end
    
    switch bg_col
        case 'white'
            col_node = [1 1 1].*0.1;
            col_line_all = [0.3 0.3 0.3]; % color of reactions in each layer
            col_edge_all = [0.9 0.9 0.9]; % color of edges not selected
            col_highlight = {[1 1 0],[0 1 0]};
            col_frame = 'k';
        case 'black'
            col_node = [1 1 1].*0.9;
            col_line_all = [0.8 0.8 0.8];
            col_edge_all = [0.3 0.3 0.3];
            col_highlight = {[1 0 0],[0 1 0]};
            col_frame = 'w';
    end
    
    % set default node color
    node_col_list = {col_node,[],col_node,col_node};
    node_list = {'met','flux','pro','rna'};
    
    % set default edge color
    col_edge_list = {lib.my_colors('matlab_orange',1,false),...
        lib.my_colors('matlab_yellow',1,false),...
        lib.my_colors('matlab_blue',1,false),...
        lib.my_colors('matlab_lblue',1,false),...
        lib.my_colors('matlab_red',1,false),...
        lib.my_colors('matlab_purple',1,false)};
    edge_list = {'s-flux','p-flux','pro-flux','pro-rna','co-flux','allo-flux'};
    
    % plot
    switch type
        case 'singleLayer'
            xy_col = node_col_list{ismember(name,node_list)};
            plot_single(xy_col,col_line_all,col_frame,col_highlight,fig,varargin);
            
        case 'multiLayer'
            num_layer = length(varargin{1});
            for i=1:num_layer
                args_now = cellfun(@(x) x{i}, varargin, 'UniformOutput',false);
                xy_col = node_col_list{ismember(node_list,name{i})};
                plot_single(xy_col,col_line_all,col_frame,col_highlight,fig,args_now);
            end
            
        case 'edges'
            num_edges = length(varargin{1});
            for i=1:num_edges
                args_now = cellfun(@(x) x{i}, varargin, 'UniformOutput',false);
                col_edge = col_edge_list{ismember(edge_list,name{i})};
                plot_edge(col_edge,fig,args_now);
            end
   
        case 'edges_sub'
            num_edges = length(varargin{1});
            for i=1:num_edges
                args_now = cellfun(@(x) x{i}, varargin, 'UniformOutput',false);
                col_edge = col_edge_list{ismember(edge_list,name{i})};
                plot_edge_sub(col_edge,col_edge_all,fig,args_now);
            end
    end
    
    
end

function plot_single(col_xy,col_line_all,col_frame,col_highlight,fig,args)
    
    n = 2;% resolution of line (# of points)
    xy = args{1};
    line_rxn = args{2};
    line_rxn_all = args{3};
    frame = args{4};
    data_xy = args{5};
    
    num_plot = size(xy,1);
    threshold_fc = 1.5;% threshold to highlight lines/nodes
    
    data_xy(isnan(data_xy)) = 0;
    if isempty(data_xy)
        data_xy = zeros(num_plot,1);
    end
    
    
    figure(fig);
    hold on;
    
    for r=1:size(line_rxn_all,1)
        line(linspace(line_rxn_all(r,1),line_rxn_all(r,4),n),...
            linspace(line_rxn_all(r,2),line_rxn_all(r,5),n),...
            linspace(line_rxn_all(r,3),line_rxn_all(r,6),n),...
            'Color',col_line_all,'LineWidth',1);
    end
    
    % if flux: col_xy = []
    col_flux = [1 1 1].*0.9;
    if isempty(col_xy)
        for r=1:size(line_rxn,1)
            if data_xy(r)>threshold_fc
                line(linspace(line_rxn(r,1),line_rxn(r,4),n),...
                    linspace(line_rxn(r,2),line_rxn(r,5),n),...
                    linspace(line_rxn(r,3),line_rxn(r,6),n),...
                    'Color',col_highlight{1},'LineWidth',1+data_xy(r)*2.5);
%                 scatter3( (line_rxn(r,1)+line_rxn(r,4))/2,...
%                     (line_rxn(r,2)+line_rxn(r,5))/2,...
%                     (line_rxn(r,3)+line_rxn(r,6))/2,...
%                     'filled','MarkerFaceColor',col_highlight{1},'SizeData',1+data_xy(r)*25);
            elseif  data_xy(r)<1/threshold_fc && data_xy(r)>0
                line(linspace(line_rxn(r,1),line_rxn(r,4),n),...
                    linspace(line_rxn(r,2),line_rxn(r,5),n),...
                    linspace(line_rxn(r,3),line_rxn(r,6),n),...
                    'Color',col_highlight{2},'LineWidth',1+data_xy(r)*2.5);
%                 scatter3( (line_rxn(r,1)+line_rxn(r,4))/2,...
%                     (line_rxn(r,2)+line_rxn(r,5))/2,...
%                     (line_rxn(r,3)+line_rxn(r,6))/2,...
%                     'filled','MarkerFaceColor',col_highlight{2},'SizeData',1+data_xy(r)*25);
            else
                line(linspace(line_rxn(r,1),line_rxn(r,4),n),...
                    linspace(line_rxn(r,2),line_rxn(r,5),n),...
                    linspace(line_rxn(r,3),line_rxn(r,6),n),...
                    'Color',col_flux,'LineWidth',1+data_xy(r)*2.5);
%                 scatter3( (line_rxn(r,1)+line_rxn(r,4))/2,...
%                     (line_rxn(r,2)+line_rxn(r,5))/2,...
%                     (line_rxn(r,3)+line_rxn(r,6))/2,...
%                     'filled','MarkerFaceColor',col_flux,'SizeData',1+data_xy(r)*25);
            end
        end
    end
        
    if ~isempty(xy) && ~isempty(col_xy)
        if isempty(data_xy)
            scatter3(xy(:,1),xy(:,2),xy(:,3),'filled','MarkerFaceColor',col_xy,...
            'SizeData',10);
        else
            for i=1:num_plot
                if data_xy(i)>threshold_fc
                    scatter3(xy(i,1),xy(i,2),xy(i,3),'filled',...
                        'MarkerFaceColor',col_highlight{1},'SizeData',1+data_xy(i)*25);
                elseif data_xy(i)<1/threshold_fc && data_xy(i)>0
                    scatter3(xy(i,1),xy(i,2),xy(i,3),'filled',...
                        'MarkerFaceColor',col_highlight{2},'SizeData',1+data_xy(i)*25);
                else
                    scatter3(xy(i,1),xy(i,2),xy(i,3),'filled','MarkerFaceColor',col_xy,...
                        'SizeData',1+data_xy(i)*25);
                end
            end
            
        end
        
    end
    line(frame(:,1),frame(:,2),frame(:,3),...
         'Color',col_frame,'LineWidth',2);
    
%     for f=1:4
%         ff = f+1;
%         if f==4
%            ff = 1; 
%         end
%         line(linspace(frame(f,1),frame(ff,1),n),...
%             linspace(frame(f,2),frame(ff,2),n),...
%             linspace(frame(f,3),frame(ff,3),n),...
%             'Color',col_frame,'LineWidth',2);
%     end
    
end

function plot_edge(col_edge,fig,args)
    
    n = 2;% number of linspace
    x = 7;% thickness of edges between layers
    edges = args{1};
    data_edges = args{2};
    
    if isempty(data_edges)
       data_edges = zeros(size(edges,1),1); 
    end
    
    figure(fig);
    hold on;
    
    if iscell(edges) % metabolites - flux
        for i=1:length(edges)
            edge_now = edges{i};
            if iscell(data_edges)
                data_now = data_edges{i};
            else
                data_now = data_edges(i);
            end
            if ~isempty(edge_now)
                if all(data_now==0)
                    for e=1:size(edge_now,1)
                        line(linspace(edge_now(e,1),edge_now(e,4),n),...
                            linspace(edge_now(e,2),edge_now(e,5),n),...
                            linspace(edge_now(e,3),edge_now(e,6),n),...
                            'Color',col_edge,...
                            'LineWidth',1);
                    end
                else
                    idx_nonzero = find(data_now~=0);
                    for j=1:length(idx_nonzero)
                        e = idx_nonzero(j);
                        line(linspace(edge_now(e,1),edge_now(e,4),n),...
                            linspace(edge_now(e,2),edge_now(e,5),n),...
                            linspace(edge_now(e,3),edge_now(e,6),n),...
                            'Color',col_edge,...
                            'LineWidth',x*data_now(e));
                    end
                end
            end
        end
    else % pro - rna | pro - flux
        edge_now = edges;
        data_now = data_edges;
        for e=1:size(edge_now,1)
            line(linspace(edge_now(e,1),edge_now(e,4),n),...
                linspace(edge_now(e,2),edge_now(e,5),n),...
                linspace(edge_now(e,3),edge_now(e,6),n),...
                'Color',col_edge,...
                'LineWidth',x*data_now(e));
%             1+data_now(e)*x
        end       
    end
    
end

function plot_edge_sub(edge_col,edge_col_all,fig,args)
    
    n = 2;% number of linspace
    x = 7;% thickness of edges between layers
    edges = args{1};
    data_edges = args{2};
    
    figure(fig);
    hold on;
    
    if iscell(edges)
        idx_plot = find(~cellfun(@isempty,data_edges,'UniformOutput',true));
        for i=1:length(edges)
            edge_now = edges{i};
            if any(i==idx_plot)
                data_now = data_edges{i};
            else
                data_now = [];
            end
            if ~isempty(edge_now)
                if isempty(data_now) || all(data_now==0)
                    for e=1:size(edge_now,1)
                        line(linspace(edge_now(e,1),edge_now(e,4),n),...
                            linspace(edge_now(e,2),edge_now(e,5),n),...
                            linspace(edge_now(e,3),edge_now(e,6),n),...
                            'Color',edge_col_all,...
                            'LineWidth',1);
                    end
                else
                    idx_nonzero = find(data_now~=0);
                    for j=1:length(idx_nonzero)
                        e = idx_nonzero(j);
                        if isnan(data_now(e))
                            line(linspace(edge_now(e,1),edge_now(e,4),n),...
                                linspace(edge_now(e,2),edge_now(e,5),n),...
                                linspace(edge_now(e,3),edge_now(e,6),n),...
                                'Color',edge_col_all,...
                                'LineWidth',1); 
                        else
                            line(linspace(edge_now(e,1),edge_now(e,4),n),...
                                linspace(edge_now(e,2),edge_now(e,5),n),...
                                linspace(edge_now(e,3),edge_now(e,6),n),...
                                'Color',edge_col,...
                                'LineWidth',x*(data_now(e)));
                        end
                    end
                end
            end
        end
    else
        edge_now = edges;
        data_now = data_edges;
        idx_plot = find(~isnan(data_now));
        if isempty(data_edges)
            for e=1:size(edge_now,1)
                line(linspace(edge_now(e,1),edge_now(e,4),n),...
                    linspace(edge_now(e,2),edge_now(e,5),n),...
                    linspace(edge_now(e,3),edge_now(e,6),n),...
                    'Color',edge_col_all,...
                    'LineWidth',1);
            end
        else
            for e=1:size(edge_now,1)
                if any(e==idx_plot)
                    line(linspace(edge_now(e,1),edge_now(e,4),n),...
                        linspace(edge_now(e,2),edge_now(e,5),n),...
                        linspace(edge_now(e,3),edge_now(e,6),n),...
                        'Color',edge_col,...
                        'LineWidth',x*data_now(e));
                else
                    line(linspace(edge_now(e,1),edge_now(e,4),n),...
                        linspace(edge_now(e,2),edge_now(e,5),n),...
                        linspace(edge_now(e,3),edge_now(e,6),n),...
                        'Color',edge_col_all,...
                        'LineWidth',1);
                end
            end
        end
        
    end
    
end