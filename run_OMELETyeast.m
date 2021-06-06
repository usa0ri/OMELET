%%%%%%%%%%%%%%%%%%%
% make datasets by simulation of kinetic model
addpath ./make_input_OMELET

obj = simKineticModel('Messiha2013','WT');

sbml_path = './data/';
savedir = './data/simulation_yeast';
mkdir(savedir);
obj.makeWT(sbml_path,savedir);
strain_names = {'M01','M02','M03','M04'};
cv_list = [0.4,0.6,1.4,1.6];
for i=1:length(strain_names)
    obj.makeMut(strain_names{i},50,cv_list(i),savedir);
end

for i=1:length(strain_names)
    obj.plotData(strain_names{i},savedir)
end

save([savedir '/obj.mat'],'obj');

%%%%%%%%%%%%%%%%%%%
% make input data for RStan

savedir = './OMELET_rstan/input_OMELETyeast';
S_path = './make_input_OMELET/S_OMELETyeast.csv';
obj.makeRstanInput(S_path,savedir);

%%%%%%%%%%%%%%%%%%%

% parameter estimation by RStan