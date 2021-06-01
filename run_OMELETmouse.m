%%%%%%%%%%%%%%%%%%%
% make input data for RStan

cd make_input_OMELET/
S_path = './S_OMELETmouse.csv';
must_rxn = {'Pgm2', 'Tpi1', 'Ldha', 'Gpt','Pcx','Cs','Glud1'};
savedir = './input_OMELETmouse';
mkdir(savedir);
make_mouse_model(S_path,must_rxn,savedir);
cd ../

%%%%%%%%%%%%%%%%%%%

% parameter estimation by RStan

%%%%%%%%%%%%%%%%%%%
% make Figures

data_path = './data';
rstan_path = '/home/suematsu/Git/FluxAnalysis/RStan/result/rstan20210430_model59_OGTT0h4h';
model_path = './OMELET_rstan/input_OMELETmouse/model_data.mat';

obj = OMELETmouse(rstan_path,model_path,data_path);
savedir = '/home/suematsu/Git/FluxAnalysis_result/result20210531/OMELETmouse';
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

save([savedir '/obj.mat'],'obj');
