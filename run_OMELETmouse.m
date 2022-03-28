%%%%%%%%%%%%%%%%%%%
% make input data for RStan

addpath ./make_input_OMELET
S_path = './make_input_OMELET/S_OMELETmouse.csv';
savedir = './OMELET_rstan/input_OMELETmouse_tmp';
mkdir(savedir);
% choose indeces from Index in metabolome.csv, proteome.csv, transcriptome.csv
% 1: WT in the fasting state
% 2: WT 4h after oral glucose administration
% 3: ob/ob in the fasting state
% 4: ob/ob 4h after oral glucose administration
idx_smplgrp = [1 2 3 4];
make_input_OMELETmouse(S_path,savedir);

% make .stan + initf_*.R
fname = 'OMELETmouse';
make_rstan_model(savedir,fname);

%%%%%%%%%%%%%%%%%%%

% parameter estimation by RStan

%%%%%%%%%%%%%%%%%%%
% calculate contributions of regulators to changes in metabolic flux
% between conditions
% make Figures

data_path = './data';
rstan_path = './OMELET_rstan/result/result_mouse';
model_path = './OMELET_rstan/input_OMELETmouse/model_data.mat';

obj = outputOMELETmouse(rstan_path,model_path,data_path);
savedir = './result_OMELETmouse';
mkdir(savedir);

obj.prep;

obj.makeFig2(savedir);
obj.makeFig4B(savedir);
obj.makeFig4C_4D(savedir);
obj.makeFig5B_7B_S8B(savedir);
obj.makeFig5C_7A_S8A(savedir);
obj.makeFig6BCD_7C_S8C(savedir);
obj.makeFigS1(savedir);
obj.makeFigS4(savedir);
obj.makeFigS5(savedir);
obj.makeFigS6(savedir);
obj.makeFigS7(savedir);

save([savedir '/obj_outputOMELETmouse.mat'],'obj');
