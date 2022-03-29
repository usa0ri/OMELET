function loadData(obj,idx_smpl,col_now)

    model_data = obj.model_data;
    if ~isempty(idx_smpl)
        [num_g, num_smpl_g] = size(idx_smpl);
        idx_g = nan(num_g,2);
        idx_g(:,1) = arrayfun(@(x) num_smpl_g*(x-1)+1,1:num_g,'UniformOutput',true);
        idx_g(:,2) = num_smpl_g.*(1:num_g);
    else
        num_smpl = model_data.out.num_smpl;
        idx_g = model_data.out.idx_g;
        idx_smpl = 1:num_smpl;
        num_g = model_data.out.num_g;
        num_smpl_g = arrayfun(@(x) idx_g(x,2)-idx_g(x,1)+1,1:num_g,'UniformOutput',true);
    end

    grp_names_org = cellfun(@(x) strrep(x,'ob/ob','ob'),...
        model_data.out.grp_names,'UniformOutput',false);
    model_data.idx_g = idx_g;
    model_data.idx_g_org = idx_smpl;
    model_data.num_smpl_g = num_smpl_g;
    model_data.num_g = num_g;
    model_data.num_smpl = length(idx_smpl(:));
    model_data.grp_names = grp_names_org;
    model_data.col = col_now;
    
    % reorder sample index
    smpltype = model_data.D.smpl_type;
    mousetype = unique(smpltype(ismember(smpltype(:,1),'Mouse'),2:end));
    mousetype = cellfun(@(x) strrep(x,'ob/ob','ob'),mousetype,'UniformOutput',false);
    timepoints = unique(smpltype(ismember(smpltype(:,1),'Time'),2:end),'stable');
    timepoints = sort(cellfun(@str2double,timepoints,'UniformOutput',true),'ascend');
    
    idx_smpl = [];
    grp_names = {};
    for m=1:length(mousetype)
        for t=1:length(timepoints)
            idx_smpl = [idx_smpl...
                find(contains(grp_names_org,mousetype{m}) &...
                contains(grp_names_org,num2str(timepoints(t))))];
            grp_names = [grp_names, [mousetype{m} num2str(timepoints(t)) 'h']];
        end
    end
    
    model_data.idx_g_reorder = idx_smpl;
    model_data.grp_names = grp_names;
    model_data.mousetype = mousetype;
    model_data.timepoints = timepoints;
    
    cmb = nchoosek(1:num_g,2);
    idx_cmb = [];
    col_list = [];
    for t=1:(length(timepoints)-1)
        idx_wt_t1 = find(contains(grp_names,mousetype{1}) & contains(grp_names,num2str(timepoints(t))));
        idx_ob_t1 = find(contains(grp_names,mousetype{2}) & contains(grp_names,num2str(timepoints(t))));
        idx_wt_t2 = find(contains(grp_names,mousetype{1}) & contains(grp_names,num2str(timepoints(t+1))));
        idx_ob_t2 = find(contains(grp_names,mousetype{2}) & contains(grp_names,num2str(timepoints(t+1))));
        
        % WT(t1)h-ob(t1)h + WT(t2)h-ob(t2)h
        idx_cmb_now = find((cmb(:,1)==idx_wt_t1 & cmb(:,2)==idx_ob_t1) |...
            (cmb(:,2)==idx_wt_t1 & cmb(:,1)==idx_ob_t1));
        idx_cmb_now2 = find((cmb(:,1)==idx_wt_t2 & cmb(:,2)==idx_ob_t2) |...
            (cmb(:,2)==idx_wt_t2 & cmb(:,1)==idx_ob_t2));
        idx_cmb = [idx_cmb; idx_cmb_now idx_cmb_now2];
        col_list = [col_list; [0.3 0.3 0.3]; [0.6 0.6 0.6]];
        
        % WT(t1)h-WT(t2)h + ob(t1)h-ob(t2)h
        idx_cmb_now = find((cmb(:,1)==idx_wt_t1 & cmb(:,2)==idx_wt_t2) |...
            (cmb(:,2)==idx_wt_t1 & cmb(:,1)==idx_wt_t2));
        idx_cmb_now2 = find((cmb(:,1)==idx_ob_t1 & cmb(:,2)==idx_ob_t2) |...
            (cmb(:,2)==idx_ob_t1 & cmb(:,1)==idx_ob_t2));
        idx_cmb = [idx_cmb; idx_cmb_now idx_cmb_now2];
        col_list = [col_list; model_data.col(idx_wt_t1,:); model_data.col(idx_ob_t1,:)];
    end
    
    model_data.idx_cmb = idx_cmb;
    model_data.col_cmb = col_list;
    mode_data.cmb = cmb;

    obj.model_data = model_data;

end