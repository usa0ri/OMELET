function plotData(obj,strain_name,savedir)

data = load_data(obj,strain_name);

savedir_now = [savedir '/plot_' strain_name ];
mkdir(savedir_now);
is_out = hist_data(data,savedir_now);

v_wt_ref = obj.rates;
idx_tmp = ismember(obj.strain_mut_name,strain_name);
v_mut_ref = obj.rates_mut{idx_tmp};
v_wt = obj.rates_smpl;
v_mut = obj.rates_mut_smpl{idx_tmp};
% reaction idx for normalize
switch obj.model_name
    case 'vanEunen2012'
        idx_n = 1;
    case 'Messiha2013'
        [~,idx_n] = ismember('HXT',v_wt_ref(:,1));
end
hist_flux(v_wt_ref,v_mut_ref,v_wt,v_mut,is_out,idx_n,savedir_now);

end


function data = load_data(obj,strain_name)

% info
struct = obj.struct;
met_names = struct.tbl_m.varSimbio;
enz_names = struct.tbl_p.varSimbio;

% WT
tbl_out = obj.tbl_mat.tbl_out;
met_tbl = tbl_out(struct.tbl_m.idxSimbio,:);
assert(isequal(cellstr(met_tbl{:,2}),met_names));
iter = size(met_tbl,2)-4;
met_data_wt = met_tbl{:,end-iter+1:end};%% FIXME
enz_tbl = tbl_out(struct.tbl_p.idxSimbio,:);
assert(isequal(cellstr(enz_tbl{:,2}),enz_names));
enz_data_wt = enz_tbl{:,end-iter+1:end};

met_data_wt_ref = met_tbl.Value;
enz_data_wt_ref = enz_tbl.Value;

clear tbl_out met_tbl enz_tbl;

% mutant
idx_tmp = ismember(obj.strain_mut_name,strain_name);
tbl_out = obj.tbl_mat_mut{idx_tmp};
met_tbl = tbl_out.tbl_out(struct.tbl_m.idxSimbio,:);
assert(isequal(cellstr(met_tbl{:,2}),met_names));
iter_ = size(met_tbl,2)-4;
met_data_mut = met_tbl{:,end-iter_+1:end};
enz_tbl = tbl_out.tbl_out(struct.tbl_p.idxSimbio,:);
assert(isequal(cellstr(enz_tbl{:,2}),enz_names));
enz_data_mut = enz_tbl{:,end-iter_+1:end};

met_data_mut_ref = met_tbl.Value;
enz_data_mut_ref = enz_tbl.Value;


clear tbl_out met_tbl enz_tbl;

data.met{1} = met_data_wt;
data.met{2} = met_data_mut;
data.enz{1} = enz_data_wt;
data.enz{2} = enz_data_mut;
data.met_ref{1} = met_data_wt_ref;
data.met_ref{2} = met_data_mut_ref;
data.enz_ref{1} = enz_data_wt_ref;
data.enz_ref{2} = enz_data_mut_ref;
data.met_names = met_names;
data.enz_names = enz_names;

end

function is_out = hist_data(data,savedir_now)

col = {'b','r'};
lcol = {'cyan','magenta'};

% met and enz
for s = 1:2
    
    if s==1
        dat = data.met;
        dat_ref = data.met_ref;
        mu = mean([dat{1} dat{2}],2);
        varnames = data.met_names;
        fname = 'met';
        
    elseif s==2
        dat = data.enz;
        dat_ref = data.enz_ref;
        mu = mean([dat{1} dat{2}],2);
        varnames = data.enz_names;
        fname = 'enz';

    end

    num = length(varnames); 
    num_sub = ceil(sqrt(num));

    fig = figure('visible','off');
    hold on;
    for n=1:num
        subplot(num_sub,num_sub,n);
        for i=1:2
            dat_now = dat{i};
            dat_ref_now = dat_ref{i};
            if mu(n)~=0
                dat_now2 = dat_now(n,:)./mu(n); 
                dat_ref_now2 = dat_ref_now(n,:)./mu(n); 
                if i==1% common bin width
                    bin_w = range(dat_now2)/10;
                end
                hold on;
                histogram(dat_now2,...
                    'BinWidth',bin_w,...
                    'FaceColor',col{i},'FaceAlpha',0.5);
                ax = gca;
                x_seq = ax.YLim(1):(ax.YLim(2)-ax.YLim(1))/9:ax.YLim(2);
                hold on;
                line(repmat(dat_ref_now2,1,10),x_seq,...
                    'LineWidth',1,'Color',lcol{i});
            end
        end
        title(varnames{n},'Interpreter','none');
    end
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [num_sub*2 num_sub*2];
    fig.PaperPosition = [0 0 num_sub*2 num_sub*2];
    set(findobj(gcf,'Type','Axes'),'FontSize',10);
%     saveas(gcf,[savedir_now '/hist_' fname '.png' ]);
    saveas(gcf,[savedir_now '/hist_' fname '.pdf' ]);
    close all;

    % if remove outliers
    num_wt = size(dat{1},2);
    num_mut = size(dat{2},2);
    is_out = nan(num,num_wt+num_mut);
    for i=1:2
        dat_now = dat{i};
        if i==1
            is_out(:,1:num_wt) = isoutlier(dat_now','quartiles')';
        else
            is_out(:,num_wt+1:end) = isoutlier(dat_now','quartiles')';
        end
    end


    fig = figure('visible','off');
    hold on;
    for n=1:num
        subplot(num_sub,num_sub,n);
        dat1_ = dat{1};
        dat2_ = dat{2};
        dat1_ref_ = dat_ref{1};
        dat2_ref_ = dat_ref{2};
        dat1 = dat1_(n,~is_out(n,1:num_wt));
        dat2 = dat2_(n,~is_out(n,num_wt+1:end));
        mu = mean([dat1 dat2]);
        if mu~=0
            dat1 = dat1./mu;
            dat2 = dat2./mu;
            dat1_ref = dat1_ref_./mu;
            dat2_ref = dat2_ref_./mu;
            bin_w = range(dat1)/10;
            histogram(dat1,...
                'BinWidth',bin_w,...
                'FaceColor',col{1},'FaceAlpha',0.5);
            hold on;
            ax = gca;
            x_seq = ax.YLim(1):(ax.YLim(2)-ax.YLim(1))/9:ax.YLim(2);
            hold on;
            line(repmat(dat1_ref(n),1,10),x_seq,...
                'LineWidth',1,'Color',lcol{1});
            histogram(dat2,...
                'BinWidth',bin_w,...
                'FaceColor',col{2},'FaceAlpha',0.5);
            line(repmat(dat2_ref(n),1,10),x_seq,...
                'LineWidth',1,'Color',lcol{2});
        end
        title(varnames{n},'Interpreter','none');
    end
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [num_sub*2 num_sub*2];
    fig.PaperPosition = [0 0 num_sub*2 num_sub*2];
    set(findobj(gcf,'Type','Axes'),'FontSize',10);
%     saveas(gcf,[savedir_now '/hist_' fname '_rmout.png' ]);
    saveas(gcf,[savedir_now '/hist_' fname '_rmout.pdf' ]);    
    close all;

end

end

function hist_flux(v_wt_ref,v_mut_ref,v_wt,v_mut,is_out,idx_n,savedir_now)

varnames = v_wt(:,1);
num = length(varnames); 
num_sub = ceil(sqrt(num));

v_wt_data = cell2mat(v_wt(:,2:end));
v_mut_data = cell2mat(v_mut(:,2:end));
v_wt_ref_data = cell2mat(v_wt_ref(:,2));
v_mut_ref_data = cell2mat(v_mut_ref(:,2));

% normalized
v_wt_n = v_wt_data./v_wt_ref_data(idx_n);
v_wt_ref_n = v_wt_ref_data./v_wt_ref_data(idx_n);
v_mut_n = v_mut_data./v_wt_ref_data(idx_n);
v_mut_ref_n = v_mut_ref_data./v_wt_ref_data(idx_n);

col = {'b','r'};
lcol = {'cyan','magenta'};

num_wt = size(v_wt,2)-1;
num_mut = size(v_mut,2)-1;
idx_out_wt = any(is_out(:,1:num_wt),1);
idx_out_mut = any(is_out(:,num_wt+1:end),1);
    
for s=1:4
    if s==1
       v_data = {v_wt_data,v_mut_data};
       v_ref = {v_wt_ref_data, v_mut_ref_data};
       fname = 'flux';
    elseif s==2
       v_data = {v_wt_n,v_mut_n}; 
       v_ref = {v_wt_ref_n, v_mut_ref_n};
       fname = 'flux_n';
    elseif s==3 
        v_data = {v_wt_data(:,~idx_out_wt),v_mut_data(:,~idx_out_mut)};
        v_ref = {v_wt_ref_data, v_mut_ref_data};
        fname = 'flux_rmout';
    elseif s==4
        v_data = {v_wt_n(:,~idx_out_wt),v_mut_n(:,~idx_out_mut)};
        v_ref = {v_wt_ref_n, v_mut_ref_n};
        fname = 'flux_n_rmout';
    end
    
    fig = figure('visible','off');
    hold on;
    for n=1:num
        subplot(num_sub,num_sub,n);
        for i=1:2
            dat_now = v_data{i};
            dat_now2 = dat_now(n,:);
            dat_ref_now = v_ref{i};
            dat_ref_now2 = dat_ref_now(n);
            if dat_ref_now2>0 && range(dat_now2)>0
                if i==1% common bin width
                    bin_w = range(dat_now2)/10;
                end
                hold on;
                histogram(dat_now2,...
                    'BinWidth',bin_w,...
                    'FaceColor',col{i},'FaceAlpha',0.5);
                ax = gca;
                x_seq = ax.YLim(1):(ax.YLim(2)-ax.YLim(1))/9:ax.YLim(2);
                hold on;
                line(repmat(dat_ref_now2,1,10),x_seq,...
                    'LineWidth',1,'Color',lcol{i});
            end
        end
        title(varnames{n},'Interpreter','none');
    end
    fig.PaperUnits = 'inches';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [num_sub*2 num_sub*2];
    fig.PaperPosition = [0 0 num_sub*2 num_sub*2];
    set(findobj(gcf,'Type','Axes'),'FontSize',10);
%     saveas(gcf,[savedir_now '/hist_' fname '.png' ]);
    saveas(gcf,[savedir_now '/hist_' fname '.pdf' ]);
    close all;
    
end


end