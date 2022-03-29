function loadOmicsData(obj)
    
%     metabolome, proteome, and transcriptome data
    tbl_met = readtable([obj.data_path '/metabolome.csv'],...
        'ReadRowNames',false,'ReadVariableNames',false);
    tbl_pro = readtable([obj.data_path '/proteome.csv'],...
        'ReadRowNames',false,'ReadVariableNames',false);
    tbl_rna = readtable([obj.data_path '/transcriptome.csv'],...
        'ReadRowNames',false,'ReadVariableNames',false);
    
    data_met = cellfun(@(x) str2double(x), table2array(tbl_met(5:end,2:end)),'UniformOutput',true);
    data_pro = cellfun(@(x) str2double(x), table2array(tbl_pro(5:end,2:end)),'UniformOutput',true);
    data_rna = cellfun(@(x) str2double(x), table2array(tbl_rna(5:end,2:end)),'UniformOutput',true);
    
    var_met = table2array(tbl_met(5:end,1));
    var_pro = table2array(tbl_pro(5:end,1));
    var_rna = table2array(tbl_rna(5:end,1));
    
    smpltype = table2array(tbl_met(1:3,2:end))';
    smpltype_grpidx = cellfun(@(x) str2double(x),table2array(tbl_met(4,2:end)),...
        'UniformOutput',true);
    idx_smplgrp = obj.model_data.out.idx_smplgrp;
    idx_import = [];
    smpltype_grpidx_import = [];
    for g=1:length(idx_smplgrp)
       idx_import = [idx_import find(smpltype_grpidx == idx_smplgrp(g))];
       smpltype_grpidx_import = [smpltype_grpidx_import...
           repmat(idx_smplgrp(g),1,sum(smpltype_grpidx == idx_smplgrp(g)))];
    end
    num_smpl = length(idx_import);
    
    % choose variables
    [~,idx_met] = ismember(obj.model_data.out.met_names,var_met);
    [~,idx_pro] = ismember(obj.model_data.out.rxn_names,var_pro);
    [~,idx_rna] = ismember(obj.model_data.out.rxn_names,var_rna);
    data_met_out = nan(length(idx_met),num_smpl);
    data_pro_out = nan(length(idx_pro),num_smpl);
    data_rna_out = nan(length(idx_rna),num_smpl);
    data_met_out(idx_met~=0,:) = data_met(nonzeros(idx_met),idx_import);
    data_pro_out(idx_pro~=0,:) = data_pro(nonzeros(idx_pro),idx_import);
    data_rna_out(idx_rna~=0,:) = data_rna(nonzeros(idx_rna),idx_import);
    var_met_out = obj.model_data.out.met_names;
    var_pro_out = obj.model_data.out.rxn_names;
    var_rna_out = obj.model_data.out.rxn_names;
    
%     %%%%%%%%%%%%%%%%%%
% %     plasma glucose and insulin after oral glucose adminsitration
%     tbl_insulin = readtable([obj.data_path '/blood_insulin.csv'],...
%         'ReadRowNames',false,'ReadVariableNames',false); 
%     tbl_glucose = readtable([obj.data_path '/blood_glucose.csv'],...
%         'ReadRowNames',false,'ReadVariableNames',false); 
%     timepoints = table2array(tbl_insulin(1,2:end));
% %     timepoints2 = table2array(tbl_glucose(1,2:end));
%     smpltype_plasma = [repmat({'WT'},1,5) repmat({'ob'},1,5)];
%     data_insulin = table2array(tbl_insulin(2:end,2:end));
%     data_glucose = table2array(tbl_glucose(2:end,2:end));
%    
    %%%%%%%%%%%%%%%%%%
%     output
    omics_data.var_omics = {var_met_out,var_pro_out,var_rna_out};
    omics_data.data_omics = {data_met_out',data_pro_out',data_rna_out'};
    omics_data.name_omics = {'metabolite','protein','transcript'};
    omics_data.smpltype = smpltype;
    omics_data.smpltype_grpidx = smpltype_grpidx;
    
%     omics_data.data_insulin = data_insulin;
%     omics_data.data_glucose = data_glucose;
%     omics_data.smpltype_plasma = smpltype_plasma;
%     omics_data.timepoitns = timepoints;
%     
    obj.omics_data = omics_data;
    
end