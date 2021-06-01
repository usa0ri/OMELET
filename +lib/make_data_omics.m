function data_omics = make_data_omics(metabolome,proteome,transcriptome,dict_omics)

idx_var_m = nonzeros(dict_omics{1}.idx);
data_m = metabolome.data(:,idx_var_m);

idx_var_p = nonzeros(dict_omics{2}.idx);
data_p = proteome.data(:,idx_var_p);
if ~all(dict_omics{2}.idx>0)
   data_p_ = nan(size(data_p,1),length(dict_omics{2}.idx));
   data_p_(:,dict_omics{2}.idx>0) = data_p;
   data_p = data_p_;
end

idx_var_t = nonzeros(dict_omics{3}.idx);
data_t = transcriptome.data(:,idx_var_t);

data_omics{1} = data_m;
data_omics{2} = data_p;
data_omics{3} = data_t;

end