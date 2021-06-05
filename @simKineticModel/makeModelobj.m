function makeModelobj(obj,sbml_path,savedir)

model_name = obj.model_name;
strain_name = obj.strain_name;

modelobj_path = [ savedir '/modelobj/modelobj_' model_name '_' strain_name '.mat'];
% if exist(modelobj_path,'file')>0
%    obj.modelobj_path = modelobj_path; 
% else
    try
        modelObj = sbmlimport(sbml_path); 
        tol_tmp = 1e-8;
        diary([savedir '/sbiosteadystate_calc_ccc_xml.txt']);
        [ success, variant_out, modelOut, exitInfo ] = sbiosteadystate(modelObj,...
            'AbsTol',tol_tmp,'RelTol',tol_tmp,'MaxStopTime',1e+3,'MinStopTime',1);
            diary off;
    catch
        success = false; 
    end
    if success
        tbl_out = variant_parser(variant_out.Content);
        name_str = strsplit( modelObj.Name, ' ');
        model_name = name_str{1};   
        mkdir([savedir '/modelobj']);
    
        save(modelobj_path,'modelObj','variant_out','tbl_out','model_name');
    else
            error('steady state calculation failed');
    end

    obj.sbml_path = sbml_path;
    out.modelObj = modelObj;
    out.variant_out = variant_out;
    out.tbl_out = tbl_out;
    obj.modelobj = out;
% end

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