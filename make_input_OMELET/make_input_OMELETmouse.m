function [X,D,out,init] = make_mouse_model(S_path,must_rxn,savedir)
%     make inputs for rstan (omics data and stoichiometric matrix)
%     
%     S_path: path to the .cxv file containing stoichiometric matrix and
%     information on cofactors and allosteric effectors
%     e.g. S_OMELETmouse.csv
%     
%     must_rxn: user-specified independent fluxes
%     
%     savedir: directry to save the results
    
%     import the .csv file
    X = parse_tbl(S_path);

%     check whether the specified set of fluxes is independent 
    [is_ok_indflux, b_ok3_flux ] = check_indflux(X,must_rxn);
    if ~is_ok_indflux
        disp(b_ok3_flux);
       error('invalid independent flux pairs. choose combination from the above;');
    end

%     calculate kernel of the stoichiometric matrix based on independent fluxes
    X = calc_kernel(X);

    D = load_omics_data(X);


    out = make_output(X,D);

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
            fprintf(fid,[fnames{i} '\t' num2str(val) '\n']);
        end
    end
    fclose(fid);

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

    a = readtable(S_path,'ReadRowNames',true,'ReadVariableNames',true);

    % is_rxn_all: 1 included in the model
    % is_irrev_all: 1 irreversible
    % is_indflux_all: 1 independent flux
    % is_include_all: 1 used for likelihood calculation
    % is_eq_all: 1 specify the kinetic equation on stan
    is_cmplx_all = ~isnan(a{'complex',:});
    is_rxn_all = a{'rxn',:}==1;
    is_irrev_all = a{'irrev',:}==1;
    is_indflux_all = a{'indflux',:}==1;
    is_include_all = a{'include',:}==1;

    % import reaction names, metabolite names
    rxn_names_all = a.Properties.VariableNames';
    idx_eff = find(contains(a.Row,'EFFECTOR'));
    idx_num = find(contains(a.Row,'NUMBER'));

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
        met_eff_list = cell(0,2);
        for i=1:num_met_eff
            if contains(met_names_eff(i),'_')
                met_names_eff(i) = extractBefore(met_names_eff(i),'_');
            end
            idx_tmp = ismember(a.Row,met_names_eff(i));
            rxn_name_tmp = rxn_names_all(a{idx_tmp,:}==1);
            for j=1:length(rxn_name_tmp)
                met_eff_list = [met_eff_list;...
                    [met_names_eff(i),rxn_name_tmp(j)]]; 
            end
        end
    end

    % rxn_names: reaction names
    % rxn_names_indflux: independent flux names
    % rxn_names_include: reaction names used for likelihood calculation
    % rxn_names_irrev: irreversible reaction
    rxn_names = rxn_names_all(is_rxn_all);
    rxn_names_indflux = rxn_names_all(is_indflux_all);
    rxn_names_include = rxn_names_all(is_include_all);
    rxn_names_irrev = rxn_names_all(is_irrev_all&is_include_all);

    [~,idx_include] = ismember(rxn_names_include,rxn_names);
    [~,idx_indflux] = ismember(rxn_names_indflux,rxn_names);
    [~,idx_irrev] = ismember(rxn_names_irrev,rxn_names);
    is_include = ismember(rxn_names,rxn_names_include);
    is_indflux = ismember(rxn_names,rxn_names_indflux);
    is_irrev = ismember(rxn_names,rxn_names_irrev);

    % enzyme complex names
    if any(is_cmplx_all)
        rxn_names_cmplx(:,1) = rxn_names_all(is_cmplx_all);
        num_rxn_cmplx = length(rxn_names_cmplx);
        idx_cmplx_all = find(is_cmplx_all);
        for i=1:num_rxn_cmplx
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
        num_rxn_cmplx = 0;
    end

    % number of reactions
    num_r = length(rxn_names);
    num_ri = length(rxn_names_indflux);
    num_rd = num_r-num_ri;
    num_include = length(rxn_names_include);
    num_irrev = length(rxn_names_irrev);

    if contains({'Pgm2','Gpd1','Eno1'},rxn_names_include)
        num_ri_wt = num_ri-3;
    elseif contains('Eno1',rxn_names_include)
        num_ri_wt = num_ri-1;
    else
        num_ri_wt = num_ri;
    end

    % Stoichiometry matrix
    S_org = double(a{1:idx_met_row-1,:});
    S_org(isnan(S_org))=0;
    is_met = any(S_org(:,is_rxn_all),2);
    met_names = met_names_all(is_met);
    num_m = length(met_names);

    % if all the contents equals 0.5 (like only ana_TCA is considered)
    % x2
    S = S_org(is_met,is_rxn_all);
    if all(abs(nonzeros(S(:)))==0.5)
       S = S.*2; 
    end

    % end metabolites
    is_end_met = sum(S~=0,2)==1 | sum(S<0,2)==sum(S~=0,2) | sum(S>0,2)==sum(S~=0,2);
    met_names_end = met_names(is_end_met);
    S_ = S(~is_end_met,:);
    met_names_int = met_names(~is_end_met);

    % number of parameter b (elasticity for products)
    % irreversible reactions do not have b
    % Cs have 2 substrates (OAA and AcCoA) so b is elasticity for substrate2
    num_Cs = contains('Cs',rxn_names_include);
    num_b = num_include-sum(is_irrev)+num_Cs;

    % output
    idx.is_indflux = is_indflux;
    idx.is_include = is_include;
    idx.is_irrev = is_irrev;
    idx.idx_include = idx_include;
    idx.idx_indflux = idx_indflux;
    idx.idx_irrev = idx_irrev;
    idx.is_met = is_met;
    idx.is_end_met = is_end_met;

    rxn.rxn_names_all = rxn_names_all;
    rxn.rxn_names = rxn_names;
    rxn.rxn_names_indflux = rxn_names_indflux;
    rxn.rxn_names_include = rxn_names_include;
    rxn.rxn_names_irrev = rxn_names_irrev;
    rxn.rxn_names_cmplx = rxn_names_cmplx;

    num.num_r = num_r;
    num.num_ri = num_ri;
    num.num_rd = num_rd;
    num.num_rc = num_include;
    num.num_irrev = num_irrev;
    num.num_ri_wt = num_ri_wt;
    num.num_m = num_m;
    num.num_met_eff = num_met_eff;
    num.num_rxn_cmplx = num_rxn_cmplx;
    num.num_b = num_b;

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
        rxn_now = X.rxn.rxn_names(idx_now);
        if ismember(X.rxn.rxn_names_indflux{i},rxn_now)
            idx_now2 = ismember(rxn_now,X.rxn.rxn_names_indflux{i});
            idx_list(i) = idx_now(idx_now2);
        end
%         disp(num2str(i));
%         disp(X.rxn.rxn_names(idx_tmp{i}));
    end

    % if different, transform the basis by *inv(N)
    N = ker(X.idx.idx_indflux,:);
    d = det(N);
    if abs(d)<0.01
       error('unable to solve inv(N)');
    end
    ker_out = ker*inv(N);
    assert( all( eye(X.num.num_ri)-ker_out(X.idx.idx_indflux,:) <1e-3, 'all') );

    % calculate model inputs
    Si = X.sto.S_(:,X.idx.is_indflux);
    Sd = X.sto.S_(:,~X.idx.is_indflux);
    Nu = [-inv(Sd)*Si; eye(X.num.num_ri)];
    Ne = [inv(Sd); zeros(X.num.num_ri,X.num.num_rd)];

    Sp = double([Sd Si]>0);
    Sm = -double([Sd Si]<0);

    rxn_names_tmp = [X.rxn.rxn_names(~X.idx.is_indflux);X.rxn.rxn_names(X.idx.is_indflux)];
    [~,idx_calc] = ismember(X.rxn.rxn_names_include,rxn_names_tmp);
    assert(isequal(X.rxn.rxn_names(X.idx.idx_include),rxn_names_tmp(idx_calc)));

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

    X.ker = ker;

end

function D = load_omics_data(X)
%     load omics data
%     modify if you use different omics data

    [data_m2, data_p2, data_t2, is_data_m, is_data_p, is_data_t] =...
        make_data_omics(X.met.met_names,X.rxn.rxn_names);

%     data shaping
%     data_m2 -> sub, pro (log transformation)
%     data_p2 -> enz
%     data_t2 -> rna
    %%%
    % max_sub = max(sum(X.sto.S_<0,1));
    max_sub = max(sum(X.sto.S<0,1));
    % max_pro = max(sum(X.sto.S_>0,1));
    max_pro = max(sum(X.sto.S>0,1));
    num_smpl = size(data_m2,2);
    sub = zeros(X.num.num_r,max_sub,num_smpl);
    pro = zeros(X.num.num_r,max_pro,num_smpl);
    enz = zeros(X.num.num_r,num_smpl);
    rna = zeros(X.num.num_r,num_smpl);
    for i=1:X.num.num_r
       S_now = X.sto.S(:,i);

       idx_sub = find(S_now<0);

       for j=1:length(idx_sub)
           if is_data_m(idx_sub(j))>0 
               sub(i,j,:) = log(data_m2(is_data_m(idx_sub(j)),:));
           end
       end
       idx_pro = find(S_now>0);
       for j=1:length(idx_pro)
           if is_data_m(idx_pro(j))>0
               pro(i,j,:) = log(data_m2(is_data_m(idx_pro(j)),:));
           end
       end

       if is_data_p(i)
           enz(i,:) = data_p2(is_data_p(i),:);
       end 

       if is_data_t(i)
           rna(i,:) = data_t2(is_data_t(i),:);
       end 
       % convert 0 to 1
       sub(sub==0) = 1;
       pro(pro==0) = 1;
    end

    % enzyme complex and effectors
    [data_m2_eff, data_p2_eff, data_t2_eff,~,~,~] =...
            make_data_omics(X.met.met_names_eff,X.rxn.rxn_names_cmplx(:,1));

    % output
    D.sub = sub;
    D.pro = pro;
    D.enz = enz;
    D.rna = rna;

    D.met_eff = log(data_m2_eff);
    D.enz_eff = data_p2_eff;
    D.rna_eff = data_t2_eff;

    D.max_sub = max_sub;
    D.max_pro = max_pro;

    D.num_smpl = num_smpl;
    D.num_wt = 11;
    D.num_ob = 12;

end

function data = normalize_data(data,smpltype_grpidx)

    num_rxn = size(data,1);
    num_g = length(unique(smpltype_grpidx));
    for i=1:num_rxn
        for g=1:num_g
            data_now = data(i,smpltype_grpidx==g);
            data_now(data_now==0 | isnan(data_now)) = nanmean(data_now);
            data(i,smpltype_grpidx==g) = data_now;
        end
       data(i,:) = data(i,:)./nanmean(data(i,:),2);
    end
    data(data==0) = 1;

end

function [data_m2,data_p2,data_t2,is_data_m,is_data_p,is_data_t] = make_data_omics(met_names,rxn_names)
%     load omics data
%     change this function if you use different omics data
    
%%%%%%%%%%%%%%%%
%     metabolome data
    tbl_met = readtable('./data/metabolome.csv',...
        'ReadRowNames',false,'ReadVariableNames',false);
    met_names_tmp = tbl_met.Var1(5:end);
    met_names_tmp{ismember(met_names_tmp,'G3P')} = 'Glycerol 3P';
    met_names_tmp{ismember(met_names_tmp,'Fumarate')} = 'Fum';
    met_names_tmp{ismember(met_names_tmp,'Citrate')} = 'Cit';
    met_names_tmp{ismember(met_names_tmp,'Succinate')} = 'Suc';
    met_names_tmp{ismember(met_names_tmp,'Acetyl-CoA')} = 'AcCoA';
    met_names_tmp{ismember(met_names_tmp,'Malate')} = 'Mal';

    met_names_ = cellfun(@(x) char(extractBefore(x,'_')), met_names, 'UniformOutput',false);
    [~,idx_now] = ismember(met_names,met_names_tmp);
    [~,idx_now2] = ismember(met_names_,met_names_tmp);
    idx_now(idx_now2>0) = nonzeros(idx_now2);

    idx_m = idx_now;
    % idx_m(idx_m>0) = idx_m_(nonzeros(idx_m));

    data_met = cellfun(@(x) str2double(x), table2array(tbl_met(5:end,2:end)),'UniformOutput',true);
    data_m = data_met(nonzeros(idx_m),:);
    
    [~,is_data_m] = ismember(met_names,met_names_tmp(nonzeros(idx_m)));

%%%%%%%%%%%%%%%%
%     proteome and transcriptome data
    tbl_pro = readtable('./data/proteome.csv',...
        'ReadRowNames',false,'ReadVariableNames',false);
    tbl_rna = readtable('./data/transcriptome.csv',...
        'ReadRowNames',false,'ReadVariableNames',false);
    pro_names_tmp = tbl_pro.Var1(5:end);
    rna_names_tmp = tbl_rna.Var1(5:end);
    data_pro = cellfun(@(x) str2double(x), table2array(tbl_pro(5:end,2:end)),'UniformOutput',true);
    data_rna = cellfun(@(x) str2double(x), table2array(tbl_rna(5:end,2:end)),'UniformOutput',true);

%     if ~isempty(rxn_names)
        % proteome
        [~,idx_p] = ismember(rxn_names,pro_names_tmp);
        if ~strcmp(rxn_names,'Sdhb')
            idx_p(ismember(rxn_names,{'Gpt','Glud1'})) = 0;
        end
        valset_p = rxn_names;
        valset_p(idx_p==0) = cellfun(@(x) [x '_noData'],valset_p(idx_p==0),'UniformOutput',false);
        
        data_p = data_pro(nonzeros(idx_p),:);

        % transcriptome
        [~,idx_t] = ismember(rxn_names,rna_names_tmp);
        valset_t = rxn_names;
        valset_t(idx_t==0) = cellfun(@(x) [x '_noData'],valset_t(idx_t==0),'UniformOutput',false);

        data_t = data_rna(nonzeros(idx_t),:);
        
        [~,is_data_p] = ismember(rxn_names,pro_names_tmp(nonzeros(idx_p)));
        [~,is_data_t] = ismember(rxn_names,rna_names_tmp(nonzeros(idx_t)));
%     end

%%%%%%%%%%%%%%%%
%     normalization
    smpltype_grpidx = cellfun(@(x) str2double(x),table2array(tbl_met(4,2:end)),...
        'UniformOutput',true);
    data_m2 = normalize_data(data_m,smpltype_grpidx);
    data_p2 = normalize_data(data_p,smpltype_grpidx);
    data_t2 = normalize_data(data_t,smpltype_grpidx);
    

end

function out = make_output(X,D)

    out.num_r = X.num.num_r;
    out.num_ri = X.num.num_ri;
    out.num_rd = X.num.num_rd;
    out.num_rc = X.num.num_rc;
    out.num_m = X.num.num_m;
    out.num_ri_wt = X.num.num_ri_wt;
    out.num_b = X.num.num_b;

    out.num_met_eff = X.num.num_met_eff;
    out.num_rxn_cmplx = X.num.num_rxn_cmplx;

    out.num_g = 4;
    out.num_p = 2;% Gpt,Glud1

    out.S = X.sto.S;
    out.S_ = X.sto.S_;
    out.N = X.ker.N;
    out.Nu = X.ker.Nu;
    out.Ne = X.ker.Ne;
    out.Sp = X.ker.Sp;
    out.Sm = X.ker.Sm;

    out.idx_calc = X.idx_calc;
    out.idx_g = [1 11;12 23;24 35;36 47];
    out.idx_p = find(ismember(X.rxn.rxn_names_include,{'Gpt','Glud1'}));

    out.rxn_names = X.rxn.rxn_names;
    out.rxn_names_include = X.rxn.rxn_names_include;
    out.rxn_names_indflux = X.rxn.rxn_names_indflux;
    out.met_names = X.met.met_names;
    out.met_names_int = X.met.met_names_int;
    out.met_names_eff = X.met.met_names_eff;

    out.idx_indflux = X.idx.idx_indflux;
    out.is_indflux = X.idx.is_indflux;
    out.idx_include = X.idx.idx_include;
    out.is_include = X.idx.is_include;
    out.idx_irrev = X.idx.idx_irrev;
    out.is_irrev = X.idx.is_irrev;

    out.num_smpl = D.num_smpl;

    out.enz = D.enz(X.idx.idx_include,:);
    out.rna = D.rna(X.idx.idx_include,:);

    for i=1:D.max_sub
        eval(['out.sub' num2str(i) '= reshape(D.sub(X.idx.idx_include,' num2str(i) ',:),out.num_rc,out.num_smpl);']);
    end
    for i=1:D.max_pro
        eval(['out.pro' num2str(i) '= reshape(D.pro(X.idx.idx_include,' num2str(i) ',:),out.num_rc,out.num_smpl);']);
    end
    out.met_eff = D.met_eff;
    out.enz_eff = D.enz_eff;
    out.rna_eff = D.rna_eff;

end

function init = set_init(X)

    % find initial value

    Sm = X.ker.Sm;
    Sp = X.ker.Sp;
    Nu = X.ker.Nu;
    Ne = X.ker.Ne;
    c_v = 0.1;
    c_e = 0.01;

    % flux_i_names
    % keyboard;
    is_nonzero = false;
    while ~is_nonzero
    %     tmp = input('initial mu_u_wt=');
        mu_u_wt = [0.1 0.4 0.3 0.1 0.7 1 0.1]';
        mu_v_wt = X.ker.N*mu_u_wt;
    %     disp(mu_v_wt);

        dxdt_wt = (abs(Sp*mu_v_wt)+abs(Sm*mu_v_wt))./2;
        sigma_e = (c_e.*dxdt_wt).^2;
        Sigma_wt = Nu*diag(c_v.*mu_u_wt).^2*Nu' + Ne*diag(sigma_e)*Ne';
        is_nonzero = abs(det(Sigma_wt))<0.1 & all(mu_v_wt>0);
    end

    init.mu_vi_wt = mu_u_wt;
    init.mu_vi_ob = mu_u_wt.*1.5;
    init.mu_v_wt = mu_v_wt(X.idx.is_include);
    init.mu_v_ob = mu_v_wt(X.idx.is_include).*1.5;

end