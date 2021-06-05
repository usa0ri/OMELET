function rateOut(obj,strain_name)

modelobj = obj.modelobj;
switch strain_name
    case 'WT'
        tbl_mat = obj.tbl_mat;
    otherwise
        idx_tmp = ismember(obj.strain_mut_name,strain_name);
        tbl_mat = obj.tbl_mat_mut{idx_tmp};
end

tmp = [true false];
for i=1:2
    is_per_smpl = tmp(i);
    rates = calc_rates(modelobj.modelObj,tbl_mat.tbl_out,is_per_smpl,obj.model_name);
    switch strain_name
        case 'WT'
            if is_per_smpl
                obj.rates_smpl = rates;
            else
                obj.rates = rates;
            end
        otherwise
            if is_per_smpl
                obj.rates_mut_smpl{idx_tmp} = rates;
            else
                obj.rates_mut{idx_tmp} = rates;
            end
            
    end
end

end

function rates = calc_rates(modelObj,tbl_out,is_per_smpl,model_name)

%%% parse from modelObj
% flux names
vnames = get(modelObj.Reactions,'Name');
% % kinetic equations
eqs = get(modelObj.Reactions,'ReactionRate');

% replace
switch model_name
    case 'Messiha2013'
        eqs_new = cell(size(eqs));
        for i=1:length(eqs)
            eqs_new{i} = replace_var(eqs{i});
        end
end

% 
% % parameter names
% pnames = get(modelObj.Parameters,'Name');
% % parameter values
% pvals = get(modelObj.Parameters,'Value');
% 
% % compartment names
% cnames = get(modelObj.Compartments,'Name');
% % capacity of compartments
% cvals = get(modelObj.Compartments,'Capacity');
% 
% % species names
% snames = get(modelObj.Species,'Name');
% % values
% sinitvals = get(modelObj.Species,'InitialAmount');
% 
% % parse parameter/compartment values
% for p=1:length(pnames)
%    eval([ pnames{p} '=' num2str(pvals{p}) ';' ]); 
% end
% if length(cnames)>1
%     for c=1:length(cnames)
%         eval([ cnames{c} '=' num2str(cvals{c}) ';' ]); 
%     end
% else
%     eval([ cnames '=' num2str(cvals) ';']);
% end

if is_per_smpl
    iter = str2double(extractAfter(tbl_out.Properties.VariableNames{end},'data'));
    v_out = nan(length(vnames),iter);
    for s = 1:iter
        for i=1:size(tbl_out,1)
            tmp_var = replace_var(tbl_out.Name{i});
            eval([tmp_var '= tbl_out.data' num2str(s) '(' num2str(i) ');']);
        end
        for i=1:length(vnames)
            switch model_name
                case 'Messiha2013'
                    var_now = who([vnames{i} '__*']);
                    for j=1:length(var_now)
                        var_now2 = char(extractAfter(var_now{j},[vnames{i} '__']));
                        eval([var_now2 '=' num2str(eval(var_now{j})) ';' ]); 
                    end
                    v = eval(eqs_new{i});
                otherwise
                    v = eval(eqs{i});
            end
            v_out(i,s) = v;
        end
%         var_lists = cellfun(@(x) replace_var(x),tbl_out.Name,'UniformOutput',false);
%         clear(var_lists{:});
    end
else
    %%% parse from variant_out
    for i=1:size(tbl_out,1)
        tmp_var = replace_var(tbl_out.Name{i});
        eval([tmp_var '= tbl_out.Value(' num2str(i) ');']);
    end
    
    % calculate rates
    v_out = nan(length(vnames),1);
    for i=1:length(vnames)
        switch model_name
            case 'Messiha2013'
                var_now = who([vnames{i} '__*']);
                for j=1:length(var_now)
                    var_now2 = char(extractAfter(var_now{j},[vnames{i} '__']));
                    eval([var_now2 '=' num2str(eval(var_now{j})) ';' ]); 
                end
                v = eval(eqs_new{i});
            otherwise
                v = eval(eqs{i});
        end
        v_out(i) = v;
    end
end

rates = [vnames num2cell(v_out)];


end

function output = replace_var(input)

input = strrep(input,'.','__');
input = strrep(input,' ','___');

if contains(input,'[TKL___(E4P:F6P)]')
   output = ['TKLE4PF6P' char(extractAfter(input,'[TKL___(E4P:F6P)]'))]; 
elseif contains(input,'[TKL___(R5P:S7P)]')
    output = ['TKLR5PS7P' char(extractAfter(input,'[TKL___(R5P:S7P)]'))]; 
elseif contains(input,'[NADPH___oxidase]')
    output = ['NADPHoxidase' char(extractAfter(input,'[NADPH___oxidase]'))];
elseif contains(input,'[E4P___sink]')
    output = ['E4Psink' char(extractAfter(input,'[E4P___sink]'))];
elseif contains(input,'[R5P___sink]')
    output = ['R5Psink' char(extractAfter(input,'[R5P___sink]'))];
else
    output = input;
end

end
