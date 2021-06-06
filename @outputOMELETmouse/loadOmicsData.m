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
    
    col_now = [     0    90/255    1.0000;...
           77/255    196/255    1.0000;...
           1.0000         0         0;...
           1.0000    128/255    130/255];
       
    %%%%%%%%%%%%%%%%%%
%     plasma glucose and insulin after oral glucose adminsitration
    tbl_insulin = readtable([obj.data_path '/blood_insulin.csv'],...
        'ReadRowNames',false,'ReadVariableNames',false); 
    tbl_glucose = readtable([obj.data_path '/blood_glucose.csv'],...
        'ReadRowNames',false,'ReadVariableNames',false); 
    timepoints = table2array(tbl_insulin(1,2:end));
%     timepoints2 = table2array(tbl_glucose(1,2:end));
    smpltype_plasma = [repmat({'WT'},1,5) repmat({'ob'},1,5)];
    data_insulin = table2array(tbl_insulin(2:end,2:end));
    data_glucose = table2array(tbl_glucose(2:end,2:end));
   
    %%%%%%%%%%%%%%%%%%
%     output
    omics_data.var_omics = {var_met,var_pro,var_rna};
    omics_data.data_omics = {data_met',data_pro',data_rna'};
    omics_data.name_omics = {'metabolite','protein','transcript'};
    omics_data.smpltype = smpltype;
    omics_data.smpltype_grpidx = smpltype_grpidx;
    omics_data.col = col_now;
    
    omics_data.data_insulin = data_insulin;
    omics_data.data_glucose = data_glucose;
    omics_data.smpltype_plasma = smpltype_plasma;
    omics_data.timepoitns = timepoints;
    
    obj.omics_data = omics_data;
    
end