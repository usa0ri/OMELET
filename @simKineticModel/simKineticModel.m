classdef simKineticModel < handle
    % simulate ODE model by simbiology toolbox
    % analyze Rstan results
    
    properties
        model_name
        % strain_name: 'WT'
        % strain_mut_name: {'M01','M02',...}
        strain_name
        strain_mut_name
        % path to SBML file
        sbml_path
        % modelObj
        modelobj
        modelobj_mut
        % struct including index and varname
        struct
        % generated sample data
        tbl_mat
        tbl_mat_mut
        % reaction fluxes
        rates
        rates_smpl
        rates_mut
        rates_mut_smpl
    end
    
    methods
        function obj = simKineticModel(model_name,strain_name)
            % constructor
            obj.model_name = model_name;
            obj.strain_name = strain_name;
        end
       
        % steady-state simulation of SBML file
        % save obj.sbml_path obj.modelobj
        makeModelobj(obj,sbml_path,savedir);
        
        % import copasi result
        % save obj.struct
        makeIdxTbl(obj,savedir);
                
        % generate data
        % save obj.tbl_mat_path
        genData(obj,iter,cv,strain_name,savedir);
        
        % extract flux info
        rateOut(obj,strain_name);
        
        mutModelobj(obj,strain_name,cv,savedir);
       
        % make Rstan input data
        makeRstanInput(obj,S_path,noise,savedir);
        
        % WT
        function makeWT(obj,sbml_path,iter,savedir)
            obj.makeModelobj(sbml_path,savedir);
            obj.makeIdxTbl(savedir);
            cv = 0.1;
            obj.genData(iter,cv,'WT',savedir);
            obj.rateOut('WT');
        end
        
        % mutant
        function makeMut(obj,strain_name,iter,cv,savedir)
            obj.mutModelobj(strain_name,cv,savedir);
            obj.genData(iter,0.1,strain_name,savedir);
            obj.rateOut(strain_name);  
            obj.plotData(strain_name,savedir);
        end
    end

end
