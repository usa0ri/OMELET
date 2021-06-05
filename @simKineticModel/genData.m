function genData(obj,iter,cv,strain_name,savedir)

savedir = [savedir '/tbl_mat'];

model_name = obj.model_name;
modelobj = obj.modelobj;
struct = obj.struct;

switch strain_name
    case 'WT'
        variant_out = modelobj.variant_out;
        tbl_out = modelobj.tbl_out;
    otherwise
        idx_tmp = ismember(obj.strain_mut_name,strain_name);
        variant_out = obj.modelobj_mut{idx_tmp}.variant_out;
        tbl_out = obj.modelobj_mut{idx_tmp}.tbl_out;
end


tbl_mat_path = [savedir '/tbl_' strain_name '_' num2str(iter) '.mat'];

    mkdir(savedir);
    modelObj = modelobj.modelObj;

    tbl_p = struct.tbl_p;
    idx_p = tbl_p.idxSimbio;
    var_p = tbl_p.varSimbio;
    tbl_m = struct.tbl_m;
    idx_m = tbl_m.idxSimbio;
    var_m = tbl_m.varSimbio;

    iter_now = 1;
    while iter_now<iter+1
        success_tmp = false;
        abs_tol = 1e-8;
        rel_tol = 1e-8;
        while ~success_tmp
            rng('shuffle');
            pert_vars = normrnd(0,cv,[1,length(var_p)]);
            v1 = sbiovariant('v1');
            v1.Content = variant_out.Content;
            for i = 1:length(var_p)
                tmp_now = v1.Content{idx_p(i)};
                pert_now = tmp_now{4}*(1+pert_vars(i));
                v1.Content{idx_p(i)}{4} = pert_now;
            end
            try
                [ success_now, varout_tmp ] = sbiosteadystate(modelObj, v1,...
                    'AbsTol',abs_tol,'RelTol',rel_tol,'MaxStopTime',1e+3,'MinStopTime',1);
                tbl_out_tmp = variant_parser(varout_tmp.Content);
                
                fc_tmp = tbl_out.Value ./ tbl_out_tmp.Value;
                fc_min = min(abs(fc_tmp));
                fc_max = max(abs(fc_tmp));
                if fc_max>2 || fc_min<0.5
                   error('too large fold change'); 
                else
                    tbl_out = addvars(tbl_out,tbl_out_tmp.Value,...
                        'NewVariableNames',['data' num2str(iter_now) ]);
                    iter_now = iter_now+1;
                end
            catch
                success_now = false;
            end
            success_tmp = success_now;
        end

    end
    disp(['dataset of ' strain_name ' with sample number ' num2str(iter) ' is successfully generated.']);
    save(tbl_mat_path,'tbl_out');
    
    noise_list = [0.01 0.05 0.1];
    idx_var = [idx_m;idx_p];
    % idx_const = find(~(tbl_out.Value==tbl_out.data1));

    data_ref = table2array(tbl_out(idx_var,4));
    data_now = table2array(tbl_out(idx_var,5:end));
    tbl_out_org = tbl_out;
    data_noise = cell(1,length(noise_list));
    i = 1;
    for n = noise_list
        rng('shuffle');
        tbl_out = tbl_out_org;
        noise_now = normrnd(0,n,size(data_now));
        noise_now_ = noise_now .* data_ref;
        data_out = abs(data_now + noise_now_);
    %     data_out2 = abs( data_now +...
    %         randn(size(data_now)).*(data_ref.*n) );
        tbl_out(idx_var,5:end) = array2table(data_out);
        save_name = [savedir '/tbl_' strain_name '_' num2str(iter) '_noise' num2str(n)];
        tbl_out_var = array2table([ idx_var data_out ],...
            'RowNames',[var_m; var_p],...
            'VariableNames',[ 'idxSimbio', tbl_out.Properties.VariableNames(5:end)]);
        save([save_name '.mat'],'tbl_out','data_out','tbl_out_var');
        
        data_noise{i} = data_out;
        i = i+1;
    end

    out.tbl_out = tbl_out_org;
    out.data_noise = data_noise;
    switch strain_name
        case 'WT'
            obj.tbl_mat = out;
        otherwise
            obj.tbl_mat_mut{idx_tmp} = out;
    end

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