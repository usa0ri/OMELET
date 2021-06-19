function outputData(obj,savedir)
    
    mkdir(savedir);
    output_data(obj,'params',savedir);
    
    output_data(obj,'contributions',savedir);
    output_data(obj,'fluxes',savedir);
    output_data(obj,'omics',savedir);
    
    
    
    
end

function output_data(obj,tbl_name,savedir)
    
    switch tbl_name
        case 'omics'
            
            % sheet names
            sheet_list = {'metabolites','proteins','transcripts',...
                'metabolitesFoldchange','proteinsFoldchange','transcriptsFoldchange'};
            
            % column names for sheet 1
            colname1 = [arrayfun(@(x) ['WT at fasting state No.' num2str(x)],1:11,'UniformOutput',false),...
                arrayfun(@(x) ['WT after oral glucose administration No.' num2str(x)],1:12,'UniformOutput',false),...
                arrayfun(@(x) ['ob/ob at fasting state No.' num2str(x)],1:12,'UniformOutput',false),...
                arrayfun(@(x) ['ob/ob after oral glucose administration No.' num2str(x)],1:12,'UniformOutput',false)];
             % column names for sheet 2
            colname2 = repmat({'Fold change','|Fold change|>1.5','q value','q value < 0.05','Change'},1,4);
            colname_list = [repmat({colname1},1,3), repmat({colname2},1,3)];

            % omics data
            [data_list,rowname_list] = calc_fc(obj);
            
            % metabolomics, proteomics, transcriptomics
            for i=1:3
                tmp = data_list{3+i};
                pmat = obj.omics_data.pmat{i};
                out = nan(size(tmp,1),size(tmp,2)*5);
                for j=1:4
                    out(:,5*j-4) = tmp(:,j);
                    out(:,5*j-3) = tmp(:,j)>1.5 | tmp(:,j)<1/1.5;
                    out(:,5*j-2:5*j-1) = pmat(2*j-1:2*j,:)';
                    out(:,5*j) = out(:,5*j-3) & out(:,5*j-1);
                end
                out2 = num2cell(out);
                for j=1:4
                    for jj=1:size(tmp,1)
                       if out(jj,5*j-3)
                           out2(jj,5*j-3) = {'TRUE'};
                       else
                           out2(jj,5*j-3) = {'FALSE'};
                       end
                       if out(jj,5*j-1)
                           out2(jj,5*j-1) = {'TRUE'};
                       else
                           out2(jj,5*j-1) = {'FALSE'};
                       end
                       if tmp(jj,j)>1.5 && out(jj,5*j-1)
                           out2(jj,5*j) = {'Increase'};
                       elseif tmp(jj,j)<1/1.5 && out(jj,5*j-1)
                           out2(jj,5*j) = {'Decrease'};
                       else
                          out2(jj,5*j) = {'Unchanged'};
                       end
                    end
                end
                data_list{3+i} = out2;
            end
                        

        case 'fluxes'
            
            sheet_list = {'Metabolic flux','Fold change','Difference'};
            
            % data list
            % 1. metabolic flux in each condition
            % 2. fold change of metabolic flux between conditions
            % 3. difference in metabolic flux in each condition
            data_list = cell(1,length(sheet_list));
            rowname_list = cell(size(data_list));
            colname_list = cell(size(data_list));
            
            num_rc = obj.model_data.out.num_rc;
            grp_names = obj.model_data.grp_names;
            num_g = obj.model_data.num_g;
            rxn_names = obj.model_data.out.rxn_names_include;
            
            % 1. metabolic flux in each condition
            v_tmp = obj.par.v;
            [v_stat, stat_names] = calc_stat_mcmc(v_tmp);
            
            data_tmp = {};
            for i=1:num_g
                v_stat_now = reshape(v_stat(:,i,:),length(stat_names),num_rc);
                colnames_now = cellfun(@(x) [ grp_names{i} ' (' x ')'],stat_names,'UniformOutput',false);
               data_tmp = [data_tmp;colnames_now;num2cell(v_stat_now')]; 
            end
            data_list{1} = data_tmp;
            rowname_list{1} = repmat([ {''};rxn_names],num_g,1);
            colname_list{1} = stat_names;
            
            % 2. fold change of metabolic flux between conditions
            cmb_list = {'Ob0h/WT0h','Ob4h/WT4h',...
                'WT4h/WT0h','Ob4h/Ob0h'};
            stat_names = [stat_names,{'Large increase (Fold change>1.5)'}];
            idx_cmb = [2 5 1 6];
            
            data_fc = cell(num_rc,length(stat_names),length(idx_cmb));
            colname_fc = cell(length(idx_cmb),length(stat_names));
            for i=1:length(idx_cmb)
                fc_tmp = obj.fc_wtob(:,:,idx_cmb(i));
                [fc_stat, ~] = calc_stat_mcmc(fc_tmp);
                fc_mean = fc_stat(2,:);
                fc_large = cell(1,num_rc);
                for ii=1:num_rc
                    if fc_mean(ii)>1.5
                        fc_large(ii) = {'True'};
                    else
                        fc_large(ii) = {'False'};
                    end
                end
                data_fc(:,:,i) = [num2cell(fc_stat') fc_large'];
                colname_fc(i,:) = cellfun(@(x) [cmb_list{i} ' (' x ')' ],stat_names,...
                    'UniformOutput',false);
            end
            % reshape data_now
            data_tmp = {};
            for j=1:length(idx_cmb)
               data_tmp = [data_tmp; colname_fc(j,:);...
                   data_fc(:,:,j)];
            end
            
            rowname_list{2} = repmat([ {''}; obj.model_data.out.rxn_names_include],length(idx_cmb),1);
            colname_list{2} = stat_names;
            data_list{2} = data_tmp;
            
            % 3. difference in metabolic flux in each condition
            cmb_list = nchoosek(1:num_rc,2);
            v_diff = nan(size(v_tmp,1),num_g,size(cmb_list,1));
            rowname_tmp = cell(size(cmb_list,1),1);
            for i=1:size(cmb_list,1)
                rowname_tmp{i} = [rxn_names{cmb_list(i,2)} ' - ' rxn_names{cmb_list(i,1)}];
                for ii=1:num_g
                    v_diff(:,ii,i) = v_tmp(:,ii,cmb_list(i,2))-v_tmp(:,ii,cmb_list(i,1));
                end
            end
            [v_diff_stat, stat_names] = calc_stat_mcmc(v_diff);
            stat_names = [stat_names,{'95% CI includes 0 or not'}];
            
            data_tmp = {};
            colnames_all = {};
            for i=1:num_g
                v_stat_now = reshape(v_diff_stat(:,i,:),length(stat_names)-1,size(cmb_list,1));
                is_diff = (v_stat_now(4,:) > 0 & v_stat_now(5,:) > 0) |...
                    (v_stat_now(4,:) < 0 & v_stat_now(5,:) < 0);
                str_diff = cell(size(is_diff));
                for ii=1:length(is_diff)
                    if is_diff(ii)
                        str_diff(ii) = {'True'};
                    else
                        str_diff(ii) = {'False'};
                    end
                end
                colnames_now = cellfun(@(x) [ grp_names{i} ' (' x ')'],stat_names,'UniformOutput',false);
                colnames_all = [colnames_all,colnames_now];
               data_tmp = [data_tmp, num2cell(v_stat_now)', str_diff']; 
            end
            data_list{3} = data_tmp;
            rowname_list{3} = [rowname_tmp];
            colname_list{3} = colnames_all;
            
            
        case 'contributions'
            % contributions of regulators to changes in metabolic flux
            % between each pair of conditions (in total 4 pairs)
            % average: average of contributions
            sheet_list = {'WT0h_Ob0h','WT4h_Ob4h',...
                'WT0h_WT4h','Ob0h_Ob4h',...
                'average'};
            
            data_list = cell(1,length(sheet_list));
            rowname_list = cell(size(data_list));
            colname_list = cell(size(data_list));
            
            idx_cmb = [2 5 1 6];
            idx_s = {1,2,3,4,5:6,7:10,11};
            rxn_names = obj.model_data.out.rxn_names_include;
            num_rc = length(rxn_names);
            iter = size(obj.par.a,1);
            s_names = {'Enzyme(transcript)','Enzyme(unaccounted)',...
               'Substrate','Product','Cofactor','Allosteric effectors','Unaccounted'};
            stat_names = {'Median','Mean','SD',...
                'Lower limit of 95% credible interval',...
                'Upper limit of 95% credible interval',...
                'Main regulator (>0.25)'};

           met_eff_list = obj.model_data.X.met.met_eff_list;
           
           % average of contributions in
            % 1. gluconeogenesis + lactate and alanine metabolism
            % 2. pyruvate cycle
            % 3. glycogenolysis + TCA cycle
           [~,idx_r_gng] = ismember({'Gpi1','Fbp1','Gpd1','Pgam1','Eno1','Ldha','Gpt'},rxn_names);
           [~,idx_r_pyr] = ismember({'Pklr','Pcx','Pck1'},rxn_names);
           [~,idx_r_other] = ismember({'Pgm2','Cs','Sdha','Fh1','Mdh2','Glud1'},rxn_names);
           idx_list = {idx_r_gng,idx_r_pyr,idx_r_other};
           idx_s_cmb = {1,1:2,3:10};% enzymes and metabolites
            s_cmb_names = {'Transcripts->Flux','Enzyme->Flux','Metabolites->Flux'};
           data_avg = nan(length(idx_list),length(s_cmb_names),length(idx_cmb));
           
            % contribution
            for i=1:length(idx_cmb)
                data_now = cell(length(rxn_names),length(stat_names),length(s_names));
                colname_now = cell(length(s_names),length(stat_names));
                
                % average
                cont_now = obj.cont_flux.intgrp(:,:,:,idx_cmb(i)); 
                for ii=1:length(idx_s_cmb)
                    data_avg(:,ii,i) = cellfun(@(x) mean(sum(cont_now(x,idx_s_cmb{ii},:),2),'all'),...
                        idx_list,'UniformOutput',true);
                end

                for j=1:length(s_names)
                    cont_ = reshape(nansum(obj.cont_flux.intgrp(:,idx_s{j},:,idx_cmb(i)),2),...
                        num_rc,iter);
                    [cont_sum_stat, ~] = calc_stat_mcmc(cont_');
                    
                    for jj=1:length(rxn_names)
                        if cont_sum_stat(2,jj)>0.25
                            data_now(jj,6,j) = {'TRUE'};
                        else
                            data_now(jj,6,j) = {'FALSE'};
                        end
                    end
                    cont_sum_stat = reshape(cont_sum_stat',num_rc,1,5);
                    cont_sum_median = cont_sum_stat(:,1,1);
                    
                    % output contributions of individual molecules
                    if j==5% cofactors
                        cont_tmp = obj.cont_flux.intgrp(:,idx_s{j},:,idx_cmb(i));
                        num_s = size(cont_tmp,2);
                        cont_tmp_stat = nan(num_rc,num_s,5);
                        for s=1:num_s
                            cont_tmp_stat(:,s,:) = calc_stat_mcmc(...
                                reshape(cont_tmp(:,s,:),num_rc,iter)')';
                        end
                        
                        cofactor_list = cellfun(@(x)...
                            met_eff_list(ismember(met_eff_list(:,2),x),1),...
                            rxn_names,'UniformOutput',false);
                        for jj=1:length(rxn_names)
                           if jj==3% Fbp1
                               cofactor_list{jj} = '';
                           elseif jj==7% Pklr
                               cofactor_list{jj} = 'ADP';
                           elseif jj==8% Pck1
                               cofactor_list{jj} = 'GTP';
                           elseif jj==11% Pcx
                               cofactor_list{jj} = 'ATP';
                           elseif jj==16% Glud1
                               cofactor_list{jj} = {'NAD+','NADH'};
                           end
                           num_eff = sum(cont_tmp_stat(jj,:,2)>0);
                           if num_eff>1
                               for jjj=1:length(stat_names)-1
                                   stat_str = arrayfun(@(x) num2str(x,'%.2f'),...
                                       cont_tmp_stat(jj,:,jjj),'UniformOutput',false);
                                   names_stat_str = arrayfun(@(x) strjoin([cofactor_list{jj}(x) stat_str(x)],' '),...
                                       1:num_eff,'UniformOutput',false);
                                   data_now(jj,jjj,j) = {[num2str(cont_sum_median(jj),'%.2f')...
                                       ' (' strjoin(names_stat_str,' + ') ')' ]};
                                   
                               end
                           end
                        end
                        
                        
                    elseif j==6% allosteric effectors
                        cont_tmp = obj.cont_flux.intgrp(:,idx_s{j},:,idx_cmb(i));
                        num_s = size(cont_tmp,2);
                        cont_tmp_stat = nan(num_rc,num_s,5);
                        for s=1:num_s
                            cont_tmp_stat(:,s,:) = calc_stat_mcmc(...
                                reshape(cont_tmp(:,s,:),num_rc,iter)')';
                        end
                        
                        allo_list = cellfun(@(x)...
                            met_eff_list(ismember(met_eff_list(:,2),x),1),...
                            rxn_names,'UniformOutput',false);
                        for jj=1:length(rxn_names)
                           allo_list{jj} = strjoin(allo_list{jj},'_');
                           if jj==3% Fbp1
                               allo_list{jj} = {'Cit(activator)','AMP(ihnibitor)'};
                           elseif any(jj==[4 8 9 13 15])% Gpd1,Pck1,Ldha,Sdha,Mdh2
                               allo_list{jj} = '';
                           elseif jj==7% Pklr
                               allo_list{jj} = {'F1,6P(activator)','ATP(inhibitor)','Ala(inhibitor)','Phe(inhibitor)'};
                           elseif jj==11% Pcx
                               allo_list{jj} = 'AcCoA(activator)';
                           elseif jj==16% Glud1
                               allo_list{jj} = {'ADP(activator)','ATP(inhibitor)','GTP(inhibitor)','Leu(inhibitor)'};
                           end
                           
                           num_eff = sum(cont_tmp_stat(jj,:,2)>0);
                           if num_eff>1
                               for jjj=1:length(stat_names)-1
                                   stat_str = arrayfun(@(x) num2str(x,'%.2f'),...
                                       cont_tmp_stat(jj,1:num_eff,jjj),'UniformOutput',false);
                                   names_stat_str = arrayfun(@(x) strjoin([allo_list{jj}(x) stat_str{x}],' '),...
                                       1:num_eff,'UniformOutput',false);
                                   data_now(jj,jjj,j) = {[num2str(cont_sum_median(jj),'%.2f')...
                                       ' (' strjoin(names_stat_str,' + ') ')' ]};
                                   
                               end
                           end
                        end
                        
                    else
                        for jjj=1:length(stat_names)-1
                            data_now(:,jjj,j) = num2cell(cont_sum_stat(:,:,jjj));
                        end
                    end
                    
                    colname_now(j,:) = cellfun(@(x) [s_names{j} ' (' x ')'],...
                        stat_names,'UniformOutput',false);

                end

                % reshape data_now
                data_tmp = {};
                for j=1:length(s_names)
                   data_tmp = [data_tmp; colname_now(j,:);...
                       data_now(:,:,j)];
                end
                
                data_list{i} = data_tmp;
                rowname_list{i} = repmat([ {''}; rxn_names],length(s_names),1);
                colname_list{i} = stat_names;
                
            end
            
            data_tmp = [];
            colnames_tmp = {};
            for i=1:length(idx_cmb)
                for ii=1:length(idx_s_cmb)
                   data_tmp = [data_tmp sum(data_avg(:,ii,i),2)]; 
                end
               colnames_tmp = [colnames_tmp,...
                   cellfun(@(x) [x ' (' sheet_list{i} ')'], s_cmb_names,'UniformOutput',false)];
            end
            data_list{end} = data_tmp;
            rowname_list{end} = {'gluconeogenesis + lactate and alanine metabolism',...
                'pyruvate cycle','glycogenolysis + TCA cycle'}';
            colname_list{end} = colnames_tmp;
            
         
        case 'params'
            
            sheet_list = {'elasticity coefficients',...
                'turnover coefficients',...
                'independent fluxes (u)',...
                'sigma',...
                'estimated amounts of OAA',...
                'estimated amounts of Pyruvate'
                };
            num_g = obj.model_data.out.num_g;
            iter = size(obj.par.a,1);
            num_rc = obj.model_data.out.num_rc;
            rxn_names = obj.model_data.out.rxn_names_include;
            grp_names = obj.model_data.grp_names;
            cmb = obj.cont.cmb;
            idx_cmb = [2 5 1 6];
            
            % parameters
            for i=1:length(sheet_list)
               data_now = [];
               name_now = {};
               if i==1
                   elasticity_type = {'Substrate or Product','Cofactor','Allosteric effectors'};
                   for ii=1:3
                      data_now = [data_now obj.elasticity_data.data{ii}];
                      name_now_tmp = cellfun(@(x) [x ' (' elasticity_type{ii} ')'],...
                          obj.elasticity_data.names{ii},'UniformOutput',false);
                      name_now = [name_now; name_now_tmp];
                   end 
               elseif i==2
                   rxn_names = [rxn_names; 'Sdhb'];
                   for g=1:num_g
                       data_now = [data_now...
                           reshape(obj.par.r_p_all(:,g,:),iter,num_rc)...
                           obj.par.r_p_eff(:,g)];
                       name_now_tmp = cellfun(@(x) [x ' (' grp_names{g} ')'],...
                          rxn_names,'UniformOutput',false);
                       name_now = [name_now; name_now_tmp];
                   end
               elseif i==3
                   num_ri = obj.model_data.out.num_ri;
                   rxn_names_indflux = obj.model_data.out.rxn_names_indflux;
                   for g=1:num_g
                      data_now = [data_now...
                          reshape(obj.par.mu_vi(:,g,:),iter,num_ri)];
                      name_now_tmp = cellfun(@(x) [x ' (' grp_names{g} ')'],...
                          rxn_names_indflux,'UniformOutput',false);
                       name_now = [name_now; name_now_tmp];
                   end
                   
               elseif i==4
                   data_now = [data_now obj.par.sigma_n obj.par.sigma_n2 obj.par.sigma_p];
                   name_now = {'sigma^{\hat e}';'sigma^{\hat t}';'sigma^{\beta}'};
                   
               elseif i==5
                   % value in each condition + fold change
                   fc_tmp = nan(size(obj.par.a,1),length(idx_cmb));
                   cmb_names = cell(1,length(idx_cmb));
                   for ii=1:length(idx_cmb)
                      fc_tmp(:,ii) = obj.par.c_OAA_out(:,cmb(idx_cmb(ii),2)) ./...
                          obj.par.c_OAA_out(:,cmb(idx_cmb(ii),1));
                      cmb_names{ii} = [grp_names{cmb(idx_cmb(ii),2)} '/' grp_names{cmb(idx_cmb(ii),1)}];
                   end
                   data_now = [data_now obj.par.c_OAA_out fc_tmp];
                   name_now = [grp_names cmb_names]';
               
               elseif i==6
                   % value in each condition + fold change
                   fc_tmp = nan(size(obj.par.a,1),length(idx_cmb));
                   cmb_names = cell(1,length(idx_cmb));
                   for ii=1:length(idx_cmb)
                      fc_tmp(:,ii) = obj.par.c_Pyr_out(:,cmb(idx_cmb(ii),2)) ./...
                          obj.par.c_Pyr_out(:,cmb(idx_cmb(ii),1));
                      cmb_names{ii} = [grp_names{cmb(idx_cmb(ii),2)} '/' grp_names{cmb(idx_cmb(ii),1)}];
                   end
                   data_now = [data_now obj.par.c_Pyr_out fc_tmp];
                   name_now = [grp_names cmb_names]';
               end

                [data_stat, stat_names] = calc_stat_mcmc(data_now);
                data_list{i} = data_stat';
                rowname_list{i} = name_now;
                colname_list{i} = stat_names;

            end
            

    end
    
    % output data
    fname = [savedir '/data_' tbl_name '.xlsx'];
    for i=1:length(sheet_list)
        if ~isempty(sheet_list{i})
            if iscell(data_list{i})
                out_now_ = [rowname_list{i}, data_list{i}];
            else
                out_now_ = [rowname_list{i},  num2cell(data_list{i})]; 
            end
           out_now = [[{''}, colname_list{i}]; out_now_];
           writecell(out_now,fname,'Sheet',sheet_list{i}); 
        end
    end
    
    
    
end

function [stat, stat_names] = calc_stat_mcmc(tmp)
    
    stat_names = {'Median','Mean','SD',...
        'Lower limit of 95% credible interval',...
        'Upper limit of 95% credible interval'};
    % calculate the following statistics from MCMC sample
    % 'Median','Mean','SD'
    % 'Lower limit of 95% credible interval'
    % 'Upper limit of 95% credible interval'
    stat_median = median(tmp,1);
    stat_mean = mean(tmp,1);
    stat_std = std(tmp,0,1);
    stat_95lower = prctile(tmp,2.5,1);
    stat_95upper = prctile(tmp,97.5,1);
    stat = cat(1,stat_median,stat_mean,stat_std,stat_95lower,stat_95upper);
    
end

function [data_list,varname_list] = calc_fc(obj)
    
    data_list = cell(1,6);
    varname_list = cell(1,6);
    idx_g = obj.model_data.idx_g;
    
    data_list{1} = obj.omics_data.data_omics{1}';
    % normalize protein data by mean of WT0h
    data_p = obj.omics_data.data_omics{2};
    data_list{2} = (data_p ./ nanmean(data_p(idx_g(1,1):idx_g(1,2),:),1))';
    data_list{3} = obj.omics_data.data_omics{3}';
    
    varname_list{1} = obj.omics_data.var_omics{1};
    varname_list{2} = obj.omics_data.var_omics{2};
    varname_list{3} = obj.omics_data.var_omics{3};
    
    % fold change
    idx_ = [1 3;2 4;1 2;3 4];
    data_list_ = cell(1,3);
    for i=1:4
        for j=1:3
            data_now = data_list{j};
            fc_now = nanmean(data_now(:,idx_g(idx_(i,2),1):idx_g(idx_(i,2),2)),2) ./...
                nanmean(data_now(:,idx_g(idx_(i,1),1):idx_g(idx_(i,1),2)),2);
            data_list_{j}{i} = fc_now;
        end
    end
    for j=1:3
       data_list{j+3} = cell2mat(data_list_{j});
    end

    varname_list{4} = obj.omics_data.var_omics{1};
    varname_list{5} = obj.omics_data.var_omics{2};
    varname_list{6} = obj.omics_data.var_omics{3};
    
end
