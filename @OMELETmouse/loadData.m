function loadData(obj,idx_smpl)

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

model_data.idx_g = idx_g;
model_data.idx_g_org = idx_smpl;
model_data.num_smpl_g = num_smpl_g;
model_data.num_g = num_g;
model_data.num_smpl = length(idx_smpl(:));
model_data.grp_names = {'WT0h','WT4h','ob0h','ob4h'};

obj.model_data = model_data;

end