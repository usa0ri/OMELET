function [X,D,out,init] = make_input_OMELETmouse(S_path,data_dir_path,idx_smplgrp,savedir)
%     make inputs for rstan (omics data and stoichiometric matrix)
%
%     import the .csv file

    % parse arguments
    if ~isfile(S_path)
        error('Error of the 1st arg: specify the path to S_*.csv');
    end
    if ~isfolder(data_dir_path)
       error('Error of the 2nd arg: specify the path to omics data'); 
    end
    if isempty(idx_smplgrp)
       error('Error of the 3rd arg: specify sample index to import'); 
    end
    if ~isfolder(savedir)
       mkdir(savedir); 
    end

    X = parse_tbl(S_path);

%     check whether the specified set of fluxes is independent 
    X = check_indflux_latest(X);

%     calculate kernel of the stoichiometric matrix based on independent fluxes
    X = calc_kernel(X);

    D = load_omics_data(X,data_dir_path,idx_smplgrp);

    out = make_output(X,D);
    
    if exist(savedir,'dir')==0
       mkdir(savedir); 
    end
    fid = fopen([savedir '/int_list.txt'],'w');
    fnames = fieldnames(out);
    for i=1:length(fnames)
        val = getfield(out,fnames{i});
        if length(val)>1
            if isnumeric(val) || islogical(val)
                writematrix(val,[savedir '/' fnames{i} '.txt']);
            else
                writecell(val,[savedir '/' fnames{i} '.txt']);
            end
        else
            if isnumeric(val)
                fprintf(fid,[fnames{i} '\t' num2str(val) '\n']);
            else
                writecell(val,[savedir '/' fnames{i} '.txt']);
            end
        end
    end
    fclose(fid);
    
    add_list = {'idx_p','idx_cmplx'};
    idx_add = find(isfield(out,add_list));
    for i=1:length(idx_add)
        if ~isempty(getfield(out,add_list{i}))
            writematrix(out.idx_p,[savedir '/' add_list{i} '.txt']);
        end
    end
    init = set_init(X);

    fnames = fieldnames(init);
    for i=1:length(fnames)
        val = getfield(init,fnames{i});
        if length(val)>1
            writematrix(val,[savedir '/' fnames{i} '.txt']);
        end
    end

    copyfile(S_path,savedir);

    save([savedir '/model_data.mat'],'X','D','out','init');

end

function X = parse_tbl(S_path)

    a_ = readtable(S_path,'ReadRowNames',true,'ReadVariableNames',true);
    % only rxn
    a = a_(:,a_{'rxn',:}==1);

    % import reaction names, metabolite names
    rxn_names_all = a.Properties.VariableNames';
    idx_eff = find(contains(a.Row,'EFFECTOR'));
    idx_num = find(contains(a.Row,'NUMBER'));
    
    number_names = a.Row(idx_num+1:end);
    assert(all(contains({'include','rxn','indflux'},number_names)));
    is_list = false(length(number_names),length(rxn_names_all));
    rxn_names_list = {};
    for i=1:length(number_names)
        is_list(i,:) = a{number_names{i},:}==1;
        rxn_names_list{i} = rxn_names_all(is_list(i,:))';
    end
    num_list = sum(is_list,2);
    % # dependent flux
    num_rd = num_list(ismember(number_names,'rxn'))-...
        num_list(ismember(number_names,'indflux'));
    % # independent fluxes to be inferred
    num_ri_wt = num_list(ismember(number_names,'indflux'))-...
        num_list(ismember(number_names,'fixed'));% since G6PC flux is fixed at 1
    % # elasticity coefficients for products (b)
    num_b = num_list(ismember(number_names,'include'))-...
        num_list(ismember(number_names,'irrev'));
    
    is_include = is_list(ismember(number_names,'include'),:);
    is_list_include = is_list(:,is_include);
    assert(all(is_list_include(ismember(number_names,'include'),:)));
    
    if isempty(idx_eff)
        idx_met_row = idx_num;
        met_names_all = a.Row(1:idx_met_row-1);
        met_eff_list = cell(0,2);
        num_met_eff = 0;
        met_names_eff = cell(0,1);
    else
        idx_met_row = idx_eff;
        met_names_all = a.Row(1:idx_met_row-1);
        met_names_eff = a.Row(idx_eff+1:idx_num-1);
        num_met_eff = length(met_names_eff);
        met_eff_list = cell(0,3);
        for i=1:num_met_eff
            idx_tmp = ismember(a.Row,met_names_eff(i));
            if contains(met_names_eff(i),'_')
                met_names_eff(i) = extractBefore(met_names_eff(i),'_');
            end
            idx_rxn = find(~isnan(a{idx_tmp,:}));
            for j=1:length(idx_rxn)
                if a{idx_tmp,idx_rxn(j)}==1
                   str_sign = '+(cofactor)';
                elseif a{idx_tmp,idx_rxn(j)}==-1
                    str_sign = '-(cofactor)';
                elseif a{idx_tmp,idx_rxn(j)}==2
                   str_sign = '+(allosteric)';
                elseif a{idx_tmp,idx_rxn(j)}==-2
                    str_sign = '-(allosteric)';
                end
                met_eff_list = [met_eff_list;...
                    [met_names_eff(i),rxn_names_all{idx_rxn(j)},str_sign]]; 
            end
        end
    end

    % enzyme complex names
    is_cmplx_all = is_list(ismember(number_names,'complex'),:);
    if any(is_cmplx_all)
        rxn_names_cmplx(:,1) = rxn_names_all(is_cmplx_all);
        num_enz_cmplx = length(rxn_names_cmplx);
        idx_cmplx_all = find(is_cmplx_all);
        for i=1:num_enz_cmplx
            is_tmp = false;
            j = idx_cmplx_all(i)-1;
            while ~is_tmp
                tmp = rxn_names_all(j);
                is_tmp = contains(tmp,rxn_names_include);
                j = j-1;
            end
            rxn_names_cmplx(i,2) = tmp;
        end
    else
        rxn_names_cmplx = cell(0,2);
        num_enz_cmplx = 0;
    end

    % Stoichiometry matrix
    S_org = double(a{1:idx_met_row-1,:});
    S_org(isnan(S_org))=0;
    S_org_ = S_org(:,is_list(ismember(number_names,'rxn'),:));
    is_met = any(S_org_,2);
    met_names = met_names_all(is_met);
    num_m = length(met_names);
    S = S_org_(is_met,:);

    % end metabolites
    is_end_met = sum(S~=0,2)==1 | sum(S<0,2)==sum(S~=0,2) | sum(S>0,2)==sum(S~=0,2);
    met_names_end = met_names(is_end_met);
    S_ = S(~is_end_met,:);
    met_names_int = met_names(~is_end_met);

    % reactions with multiple substrates/products
    is_multi_sub = find(sum(S<0,1)>1);
    disp(['Reactions with multiple substrates:' strjoin(rxn_names_all(is_multi_sub),',') newline]);
    for i=1:length(is_multi_sub)
        met_name_now = met_names(max(find(S(:,is_multi_sub(i))<0)));
       met_eff_list = [met_eff_list;...
           [met_name_now,rxn_names_all(is_multi_sub(i)),'+(substrate)']]; 
       if ~ismember(met_name_now,met_names_eff)
           met_names_eff = [met_names_eff;met_name_now];
       end
    end
    is_multi_pro = find(sum(S>0,1)>1);
    disp(['Reactions with multiple products:' strjoin(rxn_names_all(is_multi_pro),',') newline]);
    for i=1:length(is_multi_pro)
        if ~is_list(3,is_multi_pro(i))
            met_name_now = met_names(max(find(S(:,is_multi_pro(i))>0)));
            met_eff_list = [met_eff_list;...
               [met_name_now,rxn_names_all(is_multi_pro(i)),'-(product)']]; 
            if ~ismember(met_name_now,met_names_eff)
               met_names_eff = [met_names_eff;met_name_now];
            end
        end
    end
    
    % output
    for i=1:length(number_names)
       eval(['idx.is_' number_names{i} '_all=is_list(' num2str(i) ',:)'';' ]);
       eval(['idx.idx_' number_names{i} '_all=find(is_list(' num2str(i) ',:))'';' ]);
       eval(['idx.is_' number_names{i} '=is_list_include(' num2str(i) ',:)'';' ]);
       eval(['idx.idx_' number_names{i} '=find(is_list_include(' num2str(i) ',:))'';' ]);
       eval(['rxn.rxn_names_' number_names{i} '=rxn_names_list{' num2str(i) '}'';' ]);
       eval(['num.num_' number_names{i} '=num_list(' num2str(i) ');' ]);
    end
    
    idx.is_met = is_met;
    idx.is_end_met = is_end_met;
    idx.is_multi_sub = is_multi_sub;
    idx.is_multi_pro = is_multi_pro;

    num.num_rd = num_rd;
    num.num_m = num_m;
    num.num_met_eff = length(met_names_eff);
    num.num_met_eff_pairs = size(met_eff_list,1);
    num.num_b = num_b;
    num.num_ri_wt = num_ri_wt;
    num.num_mc = length(met_names_int);

    met.met_names_all = met_names_all;
    met.met_names = met_names;
    met.met_names_end = met_names_end;
    met.met_names_int = met_names_int;
    met.met_names_eff = met_names_eff;
    met.met_eff_list = met_eff_list;

    sto.S_org = S_org;
    sto.S = S;
    sto.S_ = S_;

    X.idx = idx;
    X.num = num;
    X.rxn = rxn;
    X.met = met;
    X.sto = sto;

end

function X = calc_kernel(X)

    % first trial to define independent flux
    ker = null(X.sto.S_,'r');
    assert(length(X.rxn.rxn_names_indflux)==size(ker,2));

    % check whether the defined independent fluxes are same as user definition
    idx_tmp = cell(1,size(ker,2));
    idx_list = nan(1,size(ker,2));
    for i=1:size(ker,2)
        seq_tmp = zeros(1,size(ker,2));
        seq_tmp(i) = 1;
        idx_now = find(all(ker==seq_tmp,2));
        idx_tmp{i} = idx_now;
        rxn_now = X.rxn.rxn_names_rxn(idx_now);
        if ismember(X.rxn.rxn_names_indflux{i},rxn_now)
            idx_now2 = ismember(rxn_now,X.rxn.rxn_names_indflux{i});
            idx_list(i) = idx_now(idx_now2);
        end
%         disp(num2str(i));
%         disp(X.rxn.rxn_names(idx_tmp{i}));
    end

    % if different, transform the basis by *inv(N)
    N = ker(X.idx.idx_indflux_all,:);
    d = det(N);
    if abs(d)<0.01
       error('unable to solve inv(N)');
    end
    ker_out = ker*inv(N);
    assert( all( eye(X.num.num_indflux)-ker_out(X.idx.idx_indflux_all,:) <1e-3, 'all') );

    % calculate model inputs
    Si = X.sto.S_(:,X.idx.is_indflux_all);
    Sd = X.sto.S_(:,~X.idx.is_indflux_all);
    Nu = [-inv(Sd)*Si; eye(X.num.num_indflux)];
    Ne = [inv(Sd); zeros(X.num.num_indflux,X.num.num_rd)];

    Sp = double([Sd Si]>0);
    Sm = -double([Sd Si]<0);

    rxn_names_tmp = [X.rxn.rxn_names_rxn(~X.idx.is_indflux_all);...
        X.rxn.rxn_names_rxn(X.idx.is_indflux_all)];
    [~,idx_calc] = ismember(X.rxn.rxn_names_include,rxn_names_tmp);
    assert(isequal(X.rxn.rxn_names_rxn(X.idx.idx_include_all),...
        rxn_names_tmp(idx_calc)));

    % output
    clear ker;
    ker.N = ker_out;
    ker.Si = Si;
    ker.Sd = Sd;
    ker.Nu = Nu;
    ker.Ne = Ne;
    ker.Sp = Sp;
    ker.Sm = Sm;

    X.idx_calc = idx_calc;
    
    rxn = X.rxn;
    rxn.rxn_names_Nu = [X.rxn.rxn_names_rxn(~X.idx.is_indflux_all);...
        X.rxn.rxn_names_indflux];

    X.ker = ker;
    X.rxn = rxn;

end

function D = load_omics_data(X,data_dir_path,idx_smplgrp)
%     load omics data
%     modify if you use different omics data

    [data_m2,is_data_m,cnt_nan_m,smpl_type_m] = ...
        make_data_omics([data_dir_path '/metabolome.csv'],idx_smplgrp,X.met.met_names,'metabolites');
    [data_p2,is_data_p,cnt_nan_p,smpl_type_p] = ...
        make_data_omics([data_dir_path '/proteome.csv'],idx_smplgrp,X.rxn.rxn_names_rxn,'proteins');
    [data_t2,is_data_t,cnt_nan_t,smpl_type_t] = ...
        make_data_omics([data_dir_path '/transcriptome.csv'],idx_smplgrp,X.rxn.rxn_names_rxn,'transcripts');

    assert(isequal(smpl_type_m,smpl_type_p));
    assert(isequal(smpl_type_m,smpl_type_t));
    smpl_type = smpl_type_m;

    % index of enzymes whose protein levels to be inferred from its transcript levels
    idx_p = is_data_p==0 & is_data_t>0;

%     data shaping
%     data_m2 -> sub, pro (log transformation)
%     data_p2 -> enz
%     data_t2 -> rna

    % check whether substrates and products are measured in all the
    % reaction used for likelihood calculation
    % if there are unmeasured substrates/products, infer their levels as
    % parameters
    met_est_list = cell(0,3);
    [met_est_list, max_sub] = check_measured(X,is_data_m,met_est_list,'substrate');
    [met_est_list, max_pro] = check_measured(X,is_data_m,met_est_list,'product');
    
    % prepare sub, pro, enz, rna
    num_smpl = size(smpl_type,2)-1;
    sub = ones(X.num.num_rxn,num_smpl);
    pro = ones(X.num.num_rxn,num_smpl);
    enz = zeros(X.num.num_rxn,num_smpl);
    rna = zeros(X.num.num_rxn,num_smpl);
    enz_names_est = {};
    for i=1:X.num.num_rxn
       S_now = X.sto.S(:,i);
       idx_sub = find(S_now<0);
       if is_data_m(idx_sub(1))>0 
           sub(i,:) = log(data_m2(idx_sub(1),:));
       end
       idx_pro = find(S_now>0);
       if is_data_m(idx_pro(1))>0
           pro(i,:) = log(data_m2(idx_pro(1),:));
       end

       if is_data_p(i) && is_data_t(i)
           enz(i,:) = data_p2(i,:);
           rna(i,:) = data_t2(i,:);
       elseif ~is_data_p(i) && is_data_t(i) && X.idx.is_include_all(i)
           disp('The following protein levels should be inferred from transcript levels:');
           disp([X.rxn.rxn_names_rxn{i} newline]);
           rna(i,:) = data_t2(i,:);
           enz_names_est = [enz_names_est X.rxn.rxn_names_rxn{i}];
           assert(any(i==find(idx_p)));
       end 
    end
    disp(newline);

    % enzyme complex and effectors
    % add multiple substrates/products as effectors
    if ~isempty(X.met.met_names_eff)
        [data_m_eff,is_data_m_eff,cnt_nan_m_eff,~] = ...
            make_data_omics([data_dir_path '/metabolome.csv'],idx_smplgrp,X.met.met_names_eff,'metabolite effectors');
        [met_est_list, ~, idx_eff] = check_measured(X,is_data_m_eff,met_est_list,'effector');
        met_names_est = unique(met_est_list(:,1));
        idx_est = cellfun(@(x) find(ismember(met_names_est,x)),met_est_list(:,1),'UniformOutput',true);
        num_met_est = length(met_names_est);
        
        met_eff = log(data_m_eff);
        is_met_eff_used = false(X.num.num_met_eff_pairs,1);
        patterns_tmp = {'-(cofactor)','-(product)'};
        for i=1:X.num.num_include
            idx_tmp = find(ismember(X.met.met_eff_list(:,2),X.rxn.rxn_names_include{i}));
            met_eff_tmp = X.met.met_eff_list(idx_tmp,:);
            if ~isempty(met_eff_tmp)
                for j=1:length(idx_tmp)
                   is_met_eff_used(idx_tmp(j)) = ismember(met_eff_tmp(j,1),X.met.met_names_eff(is_data_m_eff)) &...
                       ~(any(strcmp(met_eff_tmp(j,3),patterns_tmp)) && X.idx.is_irrev(i));
                end
            end
        end
        D.met_eff = met_eff;
        D.met_eff(isnan(D.met_eff)) = 0;
        D.met_eff_list_include = X.met.met_eff_list(is_met_eff_used,:);
    end
    if ~isempty(X.rxn.rxn_names_complex)
        [data_p2,is_data_p,cnt_nan_p,smpl_type_p] = ...
            make_data_omics([data_dir_path '/proteome.csv'],idx_smplgrp,X.rxn.rxn_names_rxn,'proteins');
        [data_t2,is_data_t,cnt_nan_t,smpl_type_t] = ...
            make_data_omics([data_dir_path '/transcriptome.csv'],idx_smplgrp,X.rxn.rxn_names_rxn,'transcripts');
        D.enz_eff = data_p2_eff;
        D.rna_eff = data_t2_eff;
    else
        D.enz_eff = [];
        D.rna_eff = [];
    end

    %%%% treat nan
    sub(isnan(sub)) = 0;
    pro(isnan(pro)) = 0;
    enz(isnan(enz)) = 0;
    rna(isnan(rna)) = 0;
    %%%
    
    % output
    D.sub = sub;
    D.pro = pro;
    D.enz = enz;
    D.rna = rna;

    D.max_sub = max_sub;
    D.max_pro = max_pro;

    D.num_smpl = num_smpl;
    D.num_wt = sum(ismember(smpl_type(2,:),'WT'));
    D.num_ob = sum(ismember(smpl_type(2,:),'ob/ob'));
    num_g = length(unique(smpl_type(4,2:end)));
    grp_names = cell(num_g,1);
    idx_g = nan(num_g,2);
    for g=1:num_g
        idx_tmp = find(ismember(smpl_type(4,2:end),num2str(idx_smplgrp(g))));
        idx_g(g,:) = [min(idx_tmp) max(idx_tmp)];
        grp_names{g} = [smpl_type{2,idx_tmp(end)} smpl_type{3,idx_tmp(end)} 'h'];
    end
    D.idx_g = idx_g;
    D.smpl_type = smpl_type;
    D.num_g = num_g;
    D.grp_names = grp_names;
    D.idx_smplgrp = idx_smplgrp;
    
    D.idx_p = idx_p;
    D.idx_p_ = ~idx_p;
    
    D.met_est_list = met_est_list;
    D.met_names_est = met_names_est;
    D.enz_names_est = enz_names_est;
    D.num_met_est = num_met_est;
    
end

function [met_est_list, max_met, idx_out] = check_measured(X,is_data_m,met_est_list,met_type)
    
    switch met_type
        case 'substrate'
            idx_met = arrayfun(@(x) find(X.sto.S(:,x)<0),1:X.num.num_rxn,'UniformOutput',false);
            eval_met = @(x) true;
            met_names_all = X.met.met_names;
            idx_include = X.idx.idx_include_all;
        case 'product'
            idx_met = arrayfun(@(x) find(X.sto.S(:,x)>0),1:X.num.num_rxn,'UniformOutput',false);
            eval_met = @(x) ~X.idx.is_irrev(X.idx.idx_include(x));
            met_names_all = X.met.met_names;
            idx_include = X.idx.idx_include_all;
        case 'effector'
            idx_met = cellfun(@(x)...
                find(ismember(X.met.met_names_eff,...
                X.met.met_eff_list(ismember(X.met.met_eff_list(:,2),x),1))),...
                X.rxn.rxn_names_rxn,'UniformOutput',false);
            eval_met = @(x) ~X.idx.is_irrev(X.idx.idx_include(x));
            met_names_all = X.met.met_names_eff;
            idx_include = X.idx.idx_include_all;
        otherwise
            error('Invalid met_type: met_type should be substrate, product, actitivator, inhibitor');
    end 
    
    if iscell(idx_met)
        is_measured_met = cellfun(@(x) is_data_m(x)>0,idx_met,'UniformOutput',false);
    else
        is_measured_met = arrayfun(@(x) is_data_m(x)>0,idx_met,'UniformOutput',false);
    end
    idx_met_calc = idx_met(idx_include);
    is_measured_met_calc = cellfun(@(x) x==1,is_measured_met(idx_include),'UniformOutput',false);
    max_met = max(cellfun(@(x) length(x),is_measured_met_calc,'UniformOutput',true));
    idx_out = cell(X.num.num_include,1);
    for r=1:X.num.num_include
        if iscell(idx_met_calc)
            idx_out{r} = idx_met_calc{r}(is_measured_met_calc{r});
        else
            tmp = idx_met_calc(r);
            idx_out{r} = tmp(is_measured_met_calc{r});
        end
        if any(~is_measured_met_calc{r}) && eval_met(r)
           disp(['The following ' met_type ' in reactions are not measured:']);
           rxn_names_disp = X.rxn.rxn_names_include{r};
           idx_tmp = idx_met{X.idx.idx_include_all(r)};
           met_names_disp = met_names_all{idx_tmp(~is_measured_met_calc{r})};
           disp([met_names_disp ' as ' met_type ' of ' rxn_names_disp]);
           if strcmp(met_type,'effector')
               if any(ismember(met_est_list(:,1),met_names_disp) &...
                       ismember(met_est_list(:,2),rxn_names_disp))
                    met_sign = X.met.met_eff_list{...
                        ismember(X.met.met_eff_list(:,1),met_names_disp) &...
                        ismember(X.met.met_eff_list(:,2),rxn_names_disp),3};
                    met_est_list = [met_est_list;...
                        {met_names_disp rxn_names_disp met_sign}];
               end
           else
               met_est_list = [met_est_list;...
                    {met_names_disp rxn_names_disp met_type}];
           end
        end
    end
    if size(met_est_list,1)==0
        disp(['all the ' met_type ' in reactions used for likelihood calculation are measured.']); 
    end
    %%%%%%%%%%%%%%%%%
    disp(newline);
    
end

function [data,cnt_nan] = normalize_data(data,smpltype_grpidx)

   data(data<0) = nan;
    num_rxn = size(data,1);
    grpidx = unique(smpltype_grpidx);
    num_g = length(grpidx);
    cnt_nan = zeros(num_rxn,num_g);
    for i=1:num_rxn
        for g=1:num_g
            data_now = data(i,smpltype_grpidx==grpidx(g));
            cnt_nan(i,g) = sum(data_now==0 | isnan(data_now));
            data_now(data_now==0 | isnan(data_now)) = nanmean(data_now);
            data(i,smpltype_grpidx==grpidx(g)) = data_now;
        end
       data(i,:) = data(i,:)./nanmean(data(i,:),2);
    end
    data(data==0) = 1;

end

function [data_out2,is_data_out,cnt_nan_out,smpl_type] = make_data_omics(data_path,idx_smplgrp,var_names,var_type)
    
    % FIXME: treat nan in each group??
    
    tbl_now = readtable(data_path,...
        'ReadRowNames',false,'ReadVariableNames',false);
    smpl_type = table2cell(tbl_now(1:4,:));
    var_names_tmp = tbl_now.Var1(5:end);
    smpltype_grpidx = cellfun(@(x) str2double(x),table2array(tbl_now(4,2:end)),...
        'UniformOutput',true);
    idx_import = [];
    smpltype_grpidx_import = [];
    for g=1:length(idx_smplgrp)
       idx_import = [idx_import find(smpltype_grpidx == idx_smplgrp(g))];
       smpltype_grpidx_import = [smpltype_grpidx_import...
           repmat(idx_smplgrp(g),1,sum(smpltype_grpidx == idx_smplgrp(g)))];
    end
    smpl_type = smpl_type(:,[1 1+idx_import]);

    [is_now,idx_now] = ismember(var_names,var_names_tmp);
    if ~all(is_now)
        disp(['The following ' var_type ' were not found in ' data_path ':']);
        disp([strjoin(var_names(~is_now),', ') newline]);
    end
    data_now = tbl_now(5:end,2:end).Variables;
    data_now = data_now(nonzeros(idx_now),idx_import);
    data_out = nan(length(idx_now),length(idx_import));
    data_out(is_now,:) = cellfun(@(x) str2double(x), data_now,'UniformOutput',true);
    data_out(data_out<0) = nan;
    
    %     normalization
    [data_out2,cnt_nan_out] = normalize_data(data_out,smpltype_grpidx_import);
    thres_half = arrayfun(@(x) floor(sum(smpltype_grpidx_import==x)/2),unique(smpltype_grpidx_import),...
        'UniformOutput',true);
%     thres_half = arrayfun(@(x) floor(sum(smpltype_grpidx_import==x)),unique(smpltype_grpidx_import),...
%         'UniformOutput',true);
    is_nan_out = any(cnt_nan_out>thres_half,2);
    idx_nan_out = find(is_nan_out & is_now);
%     var_names_mes = var_names(is_now);
    if ~isempty(idx_nan_out)
        disp(['The following ' var_type ' exist but too many NaN in ' data_path ':']);
        for i=1:length(idx_nan_out)
            disp([var_names{idx_nan_out(i)} ': #NaN=['...
                num2str(cnt_nan_out(idx_nan_out(i),:)) ']' ]);
        end
        disp(newline)
    end
    data_out2(idx_nan_out,:) = nan;
    is_data_out = is_now&~is_nan_out;
end

function out = make_output(X,D)

    out.num_r = X.num.num_rxn;
    out.num_ri = X.num.num_indflux;
    out.num_rd = X.num.num_rd;
    out.num_rc = X.num.num_include;
    out.num_m = X.num.num_m;
    out.num_mc = X.num.num_mc;
    out.num_ri_wt = X.num.num_ri_wt;
    out.num_b = X.num.num_b;
    out.num_smpl = D.num_smpl;

    out.num_g = D.num_g;
    out.num_p = sum(D.idx_p(X.idx.idx_include_all));
    out.num_met_eff = X.num.num_met_eff;
    out.num_met_eff_pairs = size(D.met_eff_list_include,1);
    out.num_enz_cmplx = X.num.num_complex;
    out.num_met_est = D.num_met_est;

    out.S = X.sto.S;
    out.S_ = X.sto.S_;
    out.N = X.ker.N;
    out.Nu = X.ker.Nu;
    out.Ne = X.ker.Ne;
    out.Sp = X.ker.Sp;
    out.Sm = X.ker.Sm;

    out.idx_calc = X.idx_calc;
    out.idx_g = D.idx_g;
    out.idx_p = find(D.idx_p(X.idx.idx_include_all));
    out.idx_p_ = find(D.idx_p_(X.idx.idx_include_all));
    out.idx_b = find(~ismember(X.rxn.rxn_names_include,X.rxn.rxn_names_irrev));

    out.is_indflux = X.idx.is_indflux_all;
    out.idx_indflux = X.idx.idx_indflux_all;
    out.is_include = X.idx.is_include_all;
    out.idx_include = X.idx.idx_include_all;
    out.is_irrev = X.idx.is_irrev;
    out.idx_irrev = X.idx.idx_irrev;
    out.is_cmplx = X.idx.is_complex;
    out.idx_cmplx = X.idx.idx_complex;
    out.idx_fixed = X.idx.idx_fixed_all;
    
    out.rxn_names = X.rxn.rxn_names_rxn;
    out.rxn_names_include = X.rxn.rxn_names_include;
    out.rxn_names_indflux = X.rxn.rxn_names_indflux;
    out.rxn_names_p = X.rxn.rxn_names_include(out.idx_p);
    out.rxn_names_Nu = X.rxn.rxn_names_Nu;
    out.met_names = X.met.met_names;
    out.met_names_int = X.met.met_names_int;
    out.met_names_eff = X.met.met_names_eff;
    out.met_est_list = D.met_est_list;
    out.met_names_est = D.met_names_est;
    out.met_eff_list_include = D.met_eff_list_include;

    out.enz = D.enz(X.idx.idx_include_all,:);
    out.rna = D.rna(X.idx.idx_include_all,:);
    out.sub = D.sub(X.idx.idx_include_all,:);
    out.pro = D.pro(X.idx.idx_include_all,:);
    out.met_eff = D.met_eff;

    
    out.enz_eff = D.enz_eff;
    out.rna_eff = D.rna_eff;
    
    out.enz_names_est = D.enz_names_est;
    out.grp_names = D.grp_names;
    out.idx_smplgrp = D.idx_smplgrp';

end

function init = set_init(X)

    % find initial value

    Sm = X.ker.Sm;
    Sp = X.ker.Sp;
    Nu = X.ker.Nu;
    Ne = X.ker.Ne;
    c_v = 0.1;
    c_e = 0.01;

    % calculate mu_u_wt by LP
    Aeq = Nu(X.idx.idx_fixed_all,:);
    beq = 1;
    A = -Nu;
    b = ones(size(Nu,1),1).*(-0.1);
    f = ones(1,size(Nu,2));
    lb = ones(1,size(Nu,2)).*0.1;
    ub = ones(1,size(Nu,2)).*5;
    [mu_u_wt, fmin] = linprog(f,A,b,Aeq,beq,lb,ub);
    
    % check
    mu_v_wt = X.ker.N*mu_u_wt;
    dxdt_wt = (abs(Sp*mu_v_wt)+abs(Sm*mu_v_wt))./2;
    sigma_e = (c_e.*dxdt_wt).^2;
    Sigma_wt = Nu*diag(c_v.*mu_u_wt).^2*Nu' + Ne*diag(sigma_e)*Ne';
    is_nonzero = abs(det(Sigma_wt))<0.1 & all(mu_v_wt>0);  
    assert(is_nonzero);

    init.mu_vi = mu_u_wt;
    init.mu_v = mu_v_wt(X.idx.is_include_all);

end
