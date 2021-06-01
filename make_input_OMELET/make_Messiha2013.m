function make_Messiha2013(S_path,a,savedir)

% Sinfo = load_S(a);

X = parse_tbl(S_path);

% must_rxn = {'TPS','GPD','ENO','ADH'};
must_rxn = {'TPS','GPD','acetate_branch','ADH'};
[is_ok_indflux, b_ok3_flux ] = check_indflux(X,must_rxn);
if ~is_ok_indflux
    disp(b_ok3_flux);
   error('invalid independent flux pairs. choose combination from the above;');
end

X = calc_kernel(X);

D = load_sim_data(X,a);

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

init = set_init(X,a);

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

function Sinfo = load_S(simModelObj)

load(simModelObj.modelobj_path);
S = getstoichmatrix(modelObj);
S = full(S);
load(simModelObj.struct_path);

idx_m = [1:29 59:67];
S_tmp = S(idx_m,:);
var_m_tmp = struct_idx_tbl.tbl_m_ss_copasi.Properties.RowNames(idx_m);
var_r_tmp = struct_idx_tbl.tbl_r_ss_copasi.Properties.RowNames;

Sinfo.S_full = S;
Sinfo.S = S_tmp;
Sinfo.m_names = var_m_tmp;
Sinfo.r_names = var_r_tmp;

r_tmp = matlab.lang.makeValidName(var_r_tmp);
m_tmp = matlab.lang.makeValidName(var_m_tmp);
tbl_tmp = array2table(Sinfo.S,'VariableNames',r_tmp,'RowNames',m_tmp);

Sinfo.tbl = tbl_tmp;

end

function X = parse_tbl(S_path)

a = readtable(S_path,'ReadRowNames',true,'ReadVariableNames',true);

% is_rxn_all: 1 included in the model
% is_irrev_all: 1 irreversible
% is_indflux_all: 1 independent flux
% is_include_all: 1 used for likelihood calculation
% is_eq_all: 1 specify the kinetic equation on stan
% is_cmplx_all = ~isnan(a{'complex',:});
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
    met_names_eff = cell(0,2);
    num_met_eff = 0;
else
    idx_met_row = idx_eff;
    met_names_all = a.Row(1:idx_met_row-1);
    met_names_eff(:,1) = a.Row(idx_eff+1:idx_num-1);
    num_met_eff = length(met_names_eff);
    for i=1:num_met_eff
        idx_tmp = ismember(a.Row,met_names_eff(i,1));
        rxn_eff = rxn_names_all(~isnan(a{idx_tmp,:})&is_include_all);
        met_names_eff(i,2:length(rxn_eff)+1) = rxn_eff';
        if contains(met_names_eff(i,1),'_')
            met_names_eff(i,1) = extractBefore(met_names_eff(i,1),'_');
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

% number of reactions
num_r = length(rxn_names);
num_ri = length(rxn_names_indflux);
num_rd = num_r-num_ri;
num_include = length(rxn_names_include);
num_irrev = length(rxn_names_irrev);
%%%%% FIXME %%%%%%%%
num_ri_wt = num_ri-1;

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
num_b = num_include-sum(is_irrev);

assert(num_rd==length(met_names_int));

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
% rxn.rxn_names_cmplx = rxn_names_cmplx;

num.num_r = num_r;
num.num_ri = num_ri;
num.num_rd = num_rd;
num.num_rc = num_include;
num.num_irrev = num_irrev;
num.num_ri_wt = num_ri_wt;
num.num_m = num_m;
num.num_met_eff = num_met_eff;
% num.num_rxn_cmplx = num_rxn_cmplx;
num.num_b = num_b;

met.met_names_all = met_names_all;
met.met_names = met_names;
met.met_names_end = met_names_end;
met.met_names_int = met_names_int;
met.met_names_eff = met_names_eff;

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
keyboard;

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
    disp(num2str(i));
    disp(X.rxn.rxn_names(idx_tmp{i}));
end

% if different, transform the basis by *inv(N)
N = ker(X.idx.idx_indflux,:);
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

function D = load_sim_data(X,a)

[dict_m,dict_p,dict_m_eff] = make_dict_sim(a,X);
[data_m,data_p,data_m_eff] = make_data_sim(a,...
    dict_m,dict_p,dict_m_eff);

[~,is_data_m] = ismember(X.met.met_names,dict_m.key(dict_m.idx>0));

max_sub = max(sum(X.sto.S_<0,1));
max_pro = max(sum(X.sto.S_>0,1));
num_smpl = size(data_m,2);
sub = zeros(X.num.num_rc,max_sub,num_smpl);
pro = zeros(X.num.num_rc,max_pro,num_smpl);
enz = zeros(X.num.num_rc,num_smpl);
for i=1:X.num.num_rc
   S_now = X.sto.S(:,X.idx.idx_include(i));
   
   idx_sub = find(S_now<0);
      
   for j=1:length(idx_sub)
       if is_data_m(idx_sub(j))>0 
           sub(i,j,:) = log(data_m(is_data_m(idx_sub(j)),:));
       end
   end
   idx_pro = find(S_now>0);
   for j=1:length(idx_pro)
       if is_data_m(idx_pro(j))>0
           pro(i,j,:) = log(data_m(is_data_m(idx_pro(j)),:));
       end
   end
   enz(i,:) = data_p(i,:); 
   
   % convert 0 to 1
   sub(sub==0) = 1;
   pro(pro==0) = 1;
end

met_eff = log(data_m_eff);
met_eff(met_eff==0) = 1;

% flux information
num_g = length(a.strain_mut_name)+1;
num_rc = X.num.num_rc;
num_ri = X.num.num_ri;
rxn_names = X.rxn.rxn_names_include;
rxn_i_names = X.rxn.rxn_names_indflux;
[~,idx_ind] = ismember(rxn_i_names,rxn_names);

[~,idx_tmp] = ismember(rxn_names,a.rates(:,1));
mu_v_wt = cell2mat(a.rates(idx_tmp,2));
mu_v_wt_smpl = cell2mat(a.rates_smpl(idx_tmp,2:end));

tmp = 1;
num_smpl_g = size(mu_v_wt_smpl,2);
for g=1:num_g
   idx_g(g,1) = tmp + num_smpl_g*(g-1) ;
   idx_g(g,2) = num_smpl_g*g;
end

mu_v_g = nan(num_rc,num_g-1);
mu_vi_g = nan(num_ri,num_g-1);
mu_v_g_smpl = nan(num_rc,num_smpl_g,num_g-1);
for g=1:(num_g-1)
    mu_v_tmp = cell2mat(a.rates_mut{g}(idx_tmp,2));
   mu_v_g(:,g) = mu_v_tmp;
   mu_v_g_smpl(:,:,g) = cell2mat(a.rates_mut_smpl{g}(idx_tmp,2:end));
   mu_vi_g(:,g) = mu_v_tmp(idx_ind);
   assert(all(X.ker.N*mu_vi_g(:,g)-mu_v_g(:,g)<1));
end

% output
D.sub = sub;
D.pro = pro;
D.enz = enz;

D.met_eff = met_eff;

D.max_sub = max_sub;
D.max_pro = max_pro;

D.num_smpl = num_smpl;

D.num_g = num_g;
D.idx_g = idx_g;

D.flux_ref = mu_v_wt./mu_v_wt;
D.flux_smpl = mu_v_wt_smpl./mu_v_wt;
D.flux_mut_ref = mu_v_g./mu_v_wt;
D.flux_mut_smpl = mu_v_g_smpl./repmat(mu_v_wt,1,1,num_g);

end

function [dict_m,dict_p,dict_m_eff] = make_dict_sim(simModelObj,X)

tbl_out = simModelObj.tbl_mat.tbl_out;

met_names_now = X.met.met_names_int;
[~,idx_met] = ismember(cellfun(@(x) ['cell.' x],met_names_now,'UniformOutput',false),...
    tbl_out.Name);
met_names_now(idx_met==0) = cellfun(@(x) [x '_noData'],met_names_now(idx_met==0),...
    'UniformOutput',false);
dict_m.key = X.met.met_names_int;
dict_m.val = met_names_now;
dict_m.idx = idx_met;

%%%%%%%%%%%%%%%%%%%
rxn_list = {'cell.GND1','GND';'cell.GND2','GND2';...
    'cell.RKI1','RKI';...
    'cell.RPE1','RPE';...
    'cell.SOL3','SOL';...
    'cell.TAL1','TAL';...
    'cell.TKL1','TKL_E4PF6P';...
    'cell.TKL1','TKL_R5PS7P';...
    'cell.ZWF1','ZWF';...% E4P sink R5P sink
    'HXT.Vmax','HXT';...
    'cell.HXK1','HXK';'cell.HXK2','HXK2';'cell.GLK1','GLK1';...
    'TPS.Vmax','TPS';...
    'TPP.Vmax','TPP';...
    'UGP.Vmax','UGP';...
    'cell.PGI1','PGI';...
    'PGM.Vmax','PGM';...
    'cell.PFK2','PFK';...
    'cell.FBA1','FBA';...
    'cell.TPI1','TPI';...
    'GPD.Vmax','GPD';...
    'GPP.Vmax','GPP';...
    'cell.TDH1','TDH';'cell.TDH3','TDH3';...
    'cell.PGK1','PGK';...
    'cell.GPM1','GPM';...
    'cell.ENO1','ENO';'cell.ENO2','ENO2';...
    'cell.CDC19','PYK';...
    'cell.PDC1','PDC';'cell.PDC5','PDC5';'cell.PDC6','PDC6';...
    'acetate_branch.k','acetate_branch';...% acetate_branch
    'cell.ADH1','ADH'};
is_include = ismember(rxn_list(:,2),X.rxn.rxn_names_include);
val_rxn = rxn_list(is_include,:);
is_no_data = ~ismember(X.rxn.rxn_names_include,val_rxn(:,2));
val_out = cell(size(X.rxn.rxn_names_include));
val_out(~is_no_data) = val_rxn(:,1);
val_out(is_no_data) = cellfun(@(x) [x '_noData'],X.rxn.rxn_names_include(is_no_data),...
    'UniformOutput',false);

[~,idx_rxn] = ismember(val_out,tbl_out.Name);
dict_p.key = X.rxn.rxn_names_include;
dict_p.val = val_out;
dict_p.idx = idx_rxn;

%%%%%%%%%%%%%%%%%%%
met_names_eff = X.met.met_names_eff(:,1);
[~,idx_met_eff] = ismember(cellfun(@(x) ['cell.' x],met_names_eff,'UniformOutput',false),...
    tbl_out.Name);
met_names_eff(idx_met_eff==0) = cellfun(@(x) [x '_noData'],met_names_eff(idx_met_eff==0),...
    'UniformOutput',false);
dict_m_eff.key = X.met.met_names_eff(:,1);
dict_m_eff.val = met_names_eff;
dict_m_eff.idx = idx_met_eff;

end

function [data_m,data_p,data_m_eff] = make_data_sim(a,...
    dict_m,dict_p,dict_m_eff)

tbl_WT = a.tbl_mat.tbl_out;
tbl_mut_mat = a.tbl_mat_mut;
num_g = length(a.strain_mut_name);
num_smpl = size(tbl_WT,2)-4;

data_m = ones(length(dict_m.val),num_smpl*(num_g+1));
data_p = ones(length(dict_p.val),num_smpl*(num_g+1));
data_m_eff = ones(length(dict_m_eff.val),num_smpl*(num_g+1));

data_m_tmp = tbl_WT{nonzeros(dict_m.idx),end-50+1:end};
data_p_tmp = tbl_WT{nonzeros(dict_p.idx),end-50+1:end};
data_m_eff_tmp = tbl_WT{nonzeros(dict_m_eff.idx),end-50+1:end};
for g=1:num_g
   tbl_now = tbl_mut_mat{g}.tbl_out;
   data_m_tmp = [data_m_tmp tbl_now{nonzeros(dict_m.idx),end-50+1:end}];
   data_p_tmp = [data_p_tmp tbl_now{nonzeros(dict_p.idx),end-50+1:end}];
   data_m_eff_tmp = [data_m_eff_tmp tbl_now{nonzeros(dict_m_eff.idx),end-50+1:end}];
end

data_m(dict_m.idx>0,:) = data_m_tmp./nanmean(data_m_tmp,2);
data_p(dict_p.idx>0,:) = data_p_tmp./nanmean(data_p_tmp,2);
data_m_eff(dict_m_eff.idx>0,:) = data_m_eff_tmp./nanmean(data_m_eff_tmp,2);

end

function out = make_output(X,D)

out.num_r = X.num.num_r;
out.num_ri = X.num.num_ri;
out.num_rd = X.num.num_rd;
out.num_rc = X.num.num_rc;
out.num_m = X.num.num_m;
out.num_ri_wt = X.num.num_ri_wt;
out.num_b = X.num.num_b;

out.num_g = D.num_g;
out.idx_g = D.idx_g;

out.num_met_eff = X.num.num_met_eff;

out.S = X.sto.S;
out.S_ = X.sto.S_;
out.N = X.ker.N;
out.Nu = X.ker.Nu;
out.Ne = X.ker.Ne;
out.Sp = X.ker.Sp;
out.Sm = X.ker.Sm;

out.idx_calc = X.idx_calc;

out.rxn_names = X.rxn.rxn_names;
out.rxn_names_include = X.rxn.rxn_names_include;
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

out.enz = D.enz;
for i=1:D.max_sub
    eval(['out.sub' num2str(i) '= D.sub(:,' num2str(i) ',:);']);
end
for i=1:D.max_pro
    eval(['out.pro' num2str(i) '= D.pro(:,' num2str(i) ',:);']);
end
out.met_eff = D.met_eff;

out.flux_ref = [D.flux_ref D.flux_mut_ref];
out.flux_smpl = D.flux_smpl;
for g=1:(out.num_g-1)
   eval(['out.flux_mut_smpl' num2str(g) ' = D.flux_mut_smpl(:,:,' num2str(g) ');']); 
end


end

function init = set_init(X,a)

rxn_names = X.rxn.rxn_names_include;
rxn_i_names = X.rxn.rxn_names_indflux;
[~,idx_tmp] = ismember(rxn_names,a.rates(:,1));
idx_tmp(ismember(rxn_names,'TKL_E4PF6P')) = find(ismember(a.rates(:,1),'TKL (E4P:F6P)'));
idx_tmp(ismember(rxn_names,'TKL_R5PS7P')) = find(ismember(a.rates(:,1),'TKL (R5P:S7P)'));

mu_v_wt = cell2mat(a.rates(idx_tmp,2));
[~,idx_ind] = ismember(rxn_i_names,rxn_names);
mu_vi_wt = mu_v_wt(idx_ind);

N = X.ker.N(X.idx.idx_include,:);
assert(all(N*mu_vi_wt-mu_v_wt<1));

num_g = length(a.strain_mut_name);
mu_v_g = nan(length(rxn_names),num_g);
mu_vi_g = nan(length(rxn_i_names),num_g);
for g=1:num_g
    mu_v_tmp = cell2mat(a.rates_mut{g}(idx_tmp,2));
   mu_v_g(:,g) = mu_v_tmp;
   mu_vi_g(:,g) = mu_v_tmp(idx_ind);
   assert(all(N*mu_vi_g(:,g)-mu_v_g(:,g)<1));
end

% normalize by HXT
v_norm = mu_v_wt(ismember(rxn_names,'HXT'));

init.mu_vi_wt = mu_vi_wt./v_norm;
init.mu_v_wt = mu_v_wt./v_norm;
init.mu_vi_mut = mu_vi_g./v_norm;
init.mu_v_mut = mu_v_g./v_norm;

end

