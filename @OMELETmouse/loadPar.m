function loadPar(obj)

% model_data = obj.model_data;
data_path = obj.rstan_path;
pwd_now = pwd;
cd(data_path);
flist = dir('*.txt');
flist = flist(~cellfun('isempty', {flist.date}));
cd(pwd_now);

par_list = struct();
for i=1:length(flist)
    f_now = flist(i).name;
    par_now = char(extractBefore(f_now,'.txt'));
    val_now = importdata([data_path '/' f_now]);
    sz = size(val_now);
    if i==1
        iter = sz(1);
    end
    if isstruct(val_now)
        par_list = setfield(par_list,par_now,val_now);
    elseif sz(1)==iter
        par_list = setfield(par_list,par_now,val_now);
    else
       sz_tmp = sz(1)/iter;
       val_now2 = reshape(val_now,iter,sz_tmp,sz(2));
       par_list = setfield(par_list,par_now,val_now2);
    end
end

obj.par = par_list;
obj.par_names = fieldnames(par_list);

end