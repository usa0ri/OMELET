function X = check_indflux_latest(X)

S = X.sto.S_;
ker = null(S,'r');

num_indflux = size(ker,2);
tmp = eye(num_indflux);
rxn_names_indflux = cell(num_indflux,1);
for i=1:size(ker,2)
    is_indflux = false(size(ker,1),1);
    for j=1:size(ker,1)
        is_indflux(j) = all(tmp(i,:)==ker(j,:));
    end
    idx_indflux = find(is_indflux);
    rxn_names_indflux{i} = X.rxn.rxn_names_rxn{idx_indflux(1)};
end

X.num.num_indflux = num_indflux;
X.num.num_rd = X.num.num_rxn-num_indflux;
X.num.num_ri_wt = num_indflux-1;

X.idx.is_indflux = ismember(X.rxn.rxn_names_include,rxn_names_indflux);
[~,idx_] = ismember(rxn_names_indflux,X.rxn.rxn_names_include);
X.idx.idx_indflux = idx_;
X.idx.is_indflux_all = ismember(X.rxn.rxn_names_rxn,rxn_names_indflux);
[~,idx_] = ismember(rxn_names_indflux,X.rxn.rxn_names_rxn);
X.idx.idx_indflux_all = idx_;

X.rxn.rxn_names_indflux = rxn_names_indflux;

end