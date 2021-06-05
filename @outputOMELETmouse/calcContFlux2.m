function calcContFlux2(obj)

% {T,Eu,S,P,C,A,U}
cont_out.total = intg_cont(obj.cont.total,obj.cont_rna.total);
num_g = size(obj.cont.grp,4);
for g=1:num_g
    cont_out.grp(:,:,:,g) = intg_cont(obj.cont.grp(:,:,:,g),...
        obj.cont_rna.grp(:,:,:,g));
end
num_c = size(obj.cont.intgrp,4);
for c=1:num_c
    cont_out.intgrp(:,:,:,c) = intg_cont(obj.cont.intgrp(:,:,:,c),...
        obj.cont_rna.intgrp(:,:,:,c));
end

obj.cont_flux = cont_out;

end

function cont_out = intg_cont(cont1,cont2)

num_rc = size(cont1,1);
num_r1 = size(cont1,2);
num_r2 = size(cont2,2);
iter = size(cont1,3);

% iter x num_rc x {T,Eu,S,P,U}
cont_out = nan(num_rc,num_r1+num_r2-1,iter);

cont_out(:,1,:) = cont1(:,1,:).*cont2(:,1,:);
cont_out(:,2,:) = cont1(:,1,:).*cont2(:,2,:);
cont_out(:,3:end,:) = cont1(:,2:end,:);

assert(all(~isnan(cont_out(:))));
end

function cont_out = intg_cont2(cont1,cont2)

num_rc = size(cont1,1);
num_r1 = size(cont1,2);
num_r2 = size(cont2,2);
iter = size(cont1,3);

% iter x num_rc x {T,S,P,Eu+U}
cont_out = nan(num_rc,num_r1+num_r2-2,iter);

cont_out(:,1,:) = cont1(:,1,:).*cont2(:,1,:);
cont_out(:,2:end-1,:) = cont1(:,2:end-1,:);
cont_out(:,end,:) = cont1(:,end,:) + cont1(:,1,:).*cont2(:,end,:);

assert(all(~isnan(cont_out(:))));
end