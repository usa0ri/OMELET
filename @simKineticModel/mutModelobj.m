function mutModelobj(obj,strain_name,cv,savedir)

savedir = [savedir '/modelobj'];

model_name = obj.model_name;
mut_name = strain_name;
modelobj_wt = obj.modelobj;
modelObj = modelobj_wt.modelObj;
tbl_out = modelobj_wt.tbl_out;
struct = obj.struct;

tbl_p = struct.tbl_p;
idx_p = tbl_p.idxSimbio;
var_p = tbl_p.varSimbio;
lp = length(var_p);

if ischar(cv)
    [var_diff,tbl_diff] = load_diff_param_set(obj);
end

success_tmp = false;
tol_tmp = 1e-9*10;
abs_tol = tol_tmp;
rel_tol = tol_tmp;

count = 1;
while ~success_tmp
    rng('shuffle');
    if count>100
        abs_tol = abs_tol*10;% default:1e-8
        rel_tol = rel_tol*10;
        count = 1;
    end
    if ~ischar(cv)% add noise to each parameters
        %%%%
        int_tmp = floor(cv);
        pert_vars = int_tmp + normrnd(0,cv-int_tmp,[1 lp]);
        %%%%
        v1 = sbiovariant('v1');
        v1.Content = modelobj_wt.variant_out.Content;
        for i = 1:lp
            tmp_now = v1.Content{idx_p(i)};
            pert_now = tmp_now{4}*(2^pert_vars(i));
            v1.Content{idx_p(i)}{4} = pert_now;
        end
    else% use mutants in paper
        v1 = var_diff;
    end
    
    try
        [ success_now, varout_tmp ] = sbiosteadystate(modelObj, v1,...
            'AbsTol',abs_tol,'RelTol',rel_tol,'MaxStopTime',1e+3,'MinStopTime',1);
        tbl_out_tmp = variant_parser(varout_tmp.Content);

        fc_tmp = tbl_out_tmp.Value ./ modelobj_wt.tbl_out.Value;
        fc_min = min(abs(fc_tmp));
        fc_max = max(abs(fc_tmp));
        if fc_max>4 || fc_min<0.25
           error('too large fold change'); 
        end

    catch
        success_now = false;
        count = count + 1;
    end
    success_tmp = success_now;
end
disp([ mut_name ' is generated. (abs_tol=' num2str(abs_tol)...
    ',rel_tol=' num2str(rel_tol) ')' ]);
variant_out = varout_tmp;
tbl_out = tbl_out_tmp;
modelobj_path = [ savedir '/modelobj_' model_name '_' mut_name '.mat'];
save(modelobj_path,...
    'variant_out','tbl_out','model_name', 'mut_name');

n_tmp = length(obj.strain_mut_name);
obj.strain_mut_name{n_tmp+1} = mut_name;
out.tbl_out = tbl_out;
out.variant_out = variant_out;
obj.modelobj_mut{n_tmp+1} = out;
    
end

function [tbl_out] = variant_parser(variant_content)
    l = length(variant_content);
    tbl_out = table('Size',[l,4],'VariableTypes',{'string','string','string','double'},...
        'VariableNames',{'Type','Name','Property','Value'});
    for i = 1:l
        tbl_now = cell2table(variant_content{i});
        tbl_out(i,:) = tbl_now;
    end
end

function [var_diff,tbl_diff] = load_diff_param_set(obj)

model_name = obj.model_name;
mut_name = obj.strain_name0;
modelobj0_path = obj.modelobj_path;
struct_path = obj.struct_path;

load(modelobj0_path);
load(struct_path);

v_wt = variant_out;
tbl_wt = tbl_out;

switch model_name
    case 'Messiha2013'
        name_diff = {'o2stress'};
        
        diff_list = {'TDH.kcat_TDH1',19.1200/4;...
            'TDH.kcat_TDH3',18.1620/4;...
            '[E4P sink].k',1/100;...
            '[R5P sink].k',1/100;...
            '[NADPH oxidase].k',1*100};    
    case 'vanEunen2012'
        name_diff = {'nonstarvedD01',...
            'NstarvedD01',...
            'nonstarvedD035',...
            'NstarvedD035',...
            'glucUpshift'};
        diff_list = {'Glci',0.10,0.10,0.10,0.10,0.20;...%2
            'G6P',3.80,4.22,5.38,4.45,3.80;...%3
            'F6P',0.74,0.80,1.01,0.78,0.74;...%4
            'F16P',11.80,14.63,27.03,16.39,11.80;...%5
            'TRIO',1.00,1.00,1.00,1.00,1.00;...%6
            'BPG',0.00001,0.00001,0.00001,0.00001,0.00001;...%7
            'P3G',0.69,0.96,1.09,1.00,0.69;...%10
            'P2G',0.09,0.13,0.15,0.13,0.09;...%11
            'PEP',0.10,0.12,0.11,0.12,0.10;...%12
            'PYR',2.76,3.48,5.32,3.90,2.76;...%13
            'AcAld',0.04,0.04,0.04,0.04,0.04;...%14
            'NADH',0.29,0.29,0.29,0.29,0.29;...%9
            'Vmax_glt',220,121,201,95,160;...% 72
            'Km_glt_glco',1.6,11.0,0.9,7.0,1;...% 73
            'Km_glt_glci',1.6,11.0,0.9,7.0,1;...% 74
            'Vmax_hk',285,223,258,227,213;...%76
            'Vmax_pgi',808,852,903,856,787;...%77
            'Vmax_pfk',213,165,179,93,213;...%78
            'Vmax_ald',189,153,200,161,310;...%79
            'Vmax_gapdh_f',1859,1075,1496,853,1300;...%80
            'Vmax_gapdh_r',1211,877,867,840,853;...%81
            'Km_gapdh_gap',2.48,1.15,0.39,1.41,0.21;...%82
            'Km_gapdh_nad',2.92,2.95,2.85,2.62,2.8;...%83
            'Km_gapdh_nadh',0.022,0.10,0.007,0.014,0.06;...%84
            'Km_gapdh_bpg',1.18,0.15,0.51,1.43,0.036;...%85, Keq_gapdh (does not exist)
            'Vmax_pgk_r',2670,3030,2416,1962,2512;...%86
            'Vmax_gpm',856,748,871,403,856;...%87
            'Vmax_eno',357,285,485,272,357;...%88
            'Vmax_pyk',559,636,677,480,820;...%89
            'Vmax_pdc',248,297,335,172,395;...%90
            'Vmax_adh_r',817,744,856,744,932;...%91
            'ATP',5.00,3.92,4.29,4.70,3;...
            'ADP',1.00,0.81,1.29,1.09,1;...
            'AMP',0.30,0.25,0.44,0.37,0.3;...
            'T6P',2.20,0.36,3.52,0.59,0.2;...
            'F26P',0.014,0.009,0.003,0.0014,0.014;...
            'Glco',50,50,50,50,50;...
            'ETOH',25,25,25,25,25};
        
end

assert(all(ismember(diff_list(:,1),tbl_wt.Name)));

% only specified strain for output
idx_out = ismember(name_diff,mut_name);

l_p = size(diff_list,1);
val_wt = nan(l_p,1);
idx_v1 = nan(l_p,1);

v1 = sbiovariant('v1');
v1.Content = v_wt.Content;
for j=1:l_p
    idx_now = ismember(tbl_wt.Name,diff_list(j,1));
    tmp_now = v1.Content{idx_now};
    v1.Content{idx_now}{4} = diff_list{j,find(idx_out)+1};
    idx_v1(j,1) = find(idx_now);
    val_wt(j,1) = tmp_now{4};
end
var_diff = v1;

tbl_tmp = [diff_list num2cell(val_wt) num2cell(idx_v1)];
tbl_diff = cell2table(tbl_tmp,'VariableNames',['name' name_diff 'WT' 'idx']);

end