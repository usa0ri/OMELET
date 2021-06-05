classdef outputOMELETmouse < handle
    
    properties
        rstan_path;
        omics_data;
        data_path;
        model_data;
        
        par;
        par_names;
        
        dat_cont;
        cont;
        
        dat_cont_rna;
        cont_rna;
        
        cont_flux;
        
        fc_wtob;
        
        coord_pathway;
        coord3d_pathway;
        
        fig_info;

    end
    
    methods
        function obj = outputOMELETmouse(rstan_path,model_path,data_path)
            % rstan_path: path to the directory where results of OMELET from Rstan are saved
            % model_path: path to model_data created by make_mouse_model.m
            % data_path: path to the directory where omics data (.csv) exist
            obj.rstan_path = rstan_path;
            obj.data_path = data_path;
            obj.model_data = load(model_path);
        end
        
        loadOmicsData(obj);
        function prep(obj)
            % omics_data
            loadOmicsData(obj);
           % calculate contributions of regulator to changes in metabolic flux
           % par, par_names
            loadPar(obj);
            % update model_data as indicated sample index
            loadData(obj,[]);
            % reorder reactions
            reorderRxn(obj);
            % cont
            calcCont2(obj);
            %plotContBar(obj,'metab',savedir);
            %plotContViolin(obj,'metab',savedir);
            % cont_rna
            calcContRNA2(obj);
            %plotContBar(obj,'RNA',savedir);
            %plotContViolin(obj,'RNA',savedir);
            % cont_flux
            calcContFlux2(obj);
        end
        
        % Figure 2
        % plot metabolites, proteins, and transcripts
        makeFig2(obj,savedir);
        
        % Figure 4B
        function makeFig4B(obj,savedir)
            % fold change of metabolic fluxes (ob / WT)
            calcFCWTob(obj,[savedir '/Fig4B']);
        end
        
        % Figures 4C and 4D
        function makeFig4C_4D(obj,savedir)
            % fold change of metabolic fluxes (4h / 0h)
            calcFCfrac(obj,[savedir '/Fig4C_4D']);
        end
        
        % Figures 5B, 7B, and S8B
        function makeFig5B_7B_S8B(obj,savedir)
            plotContBar(obj,'flux',[savedir '/Fig5B_7B_S8B']);
            calcL2(obj,[savedir '/Fig7B'])
        end
        
        % Figures 5C, 7A, and S8A
        function makeFig5C_7A_S8A(obj,savedir)
            plotContFlux(obj,[savedir '/Fig5C_7A_S8A']);
        end
        
        function makeFig6BCD_7C_S8C(obj,savedir)
            % 3d plot
            inter = 4.5;
            bg_col = 'black';
            makePathway3d(obj,inter,bg_col,[savedir '/Fig6BCD_7C_S8C']);
            plotPathway3d(obj,bg_col);
            saveFigUI(obj,[savedir '/Fig6BCD_7C_S8C']);
        end
        
        % Figure S1
        % plot blood glucose and insulin
        % perform PCA for metabolites, proteins, and transcripts
        makeFigS1(obj,savedir);
        
        % Figure S4
        function makeFigS4(obj,savedir)
            % plot fitting to enzyme
            plotFitEnz(obj,[savedir '/FigS4']);
        end
        
        % Figure S5
        function makeFigS5(obj,savedir)
            % compare metabolic fluxes inferred by OMELET with those in the previous studies
            cmpTracer(obj,[savedir '/FigS5']);
        end
        
        % Figure S6
        function makeFigS6(obj,savedir)
            plotContViolin(obj,'flux',[savedir '/FigS6']);
        end
        
        % Figure S7
        function makeFigS7(obj,savedir)
            % plot elasticity
            plotElasticity(obj,[savedir '/FigS7']);
        end
        
    end
end

