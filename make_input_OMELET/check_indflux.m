function [is_ok_indflux, b_ok_flux ] = check_indflux(X,must_rxn)

    S = X.sto.S_;
    idx_indflux = X.idx.idx_indflux;
    num_ri = X.num.num_ri;

    ker = null(S,'r');
    assert(length(X.rxn.rxn_names_indflux)==size(ker,2));

    % check reactions which must be included
    [~,idx_must_list] = ismember(must_rxn,X.rxn.rxn_names_indflux);
    tmp = eye(num_ri);
    tmp(idx_must_list,:) = ker(idx_indflux(idx_must_list),:);

    is_ok_inv = all(any(tmp,1));
    num=0;
    while ~is_ok_inv
        num = num+1;
        tmp = eye(num_ri);
        tmp(idx_must_list(1:end-num),:) = ker(idx_indflux(idx_must_list(1:end-num)),:);
        is_ok_inv = all(any(tmp,1));
    end
    % if no indflux is selected
    if num==length(must_rxn)
        num2=0;
        tmp(idx_must_list,:) = ker(idx_indflux(idx_must_list),:);
        is_ok_inv = all(any(tmp,1));

        while ~is_ok_inv
            num2 = num2+1;
            tmp = eye(num_ri);
            tmp(idx_must_list(num2:end),:) = ker(idx_indflux(idx_must_list(num2:end)),:);
            is_ok_inv = all(any(tmp,1));
        end
        ker2 = ker*inv(tmp);
        num_ri_now = length(idx_must_list)-num2+1;
        idx_must_rem = 1:num2-1;
        idx_list_rem = 1:num_ri;
        idx_list_rem = idx_list_rem(~ismember(idx_list_rem,idx_must_list(num2:end)));

    else
        ker2 = ker*inv(tmp);
        num_ri_now = length(idx_must_list)-num;
        idx_must_rem = num_ri_now+1:length(idx_must_list);
        idx_list_rem = 1:num_ri;
        idx_list_rem = idx_list_rem(~ismember(idx_list_rem,idx_must_list(1:end-num)));   
    end


    % list up all the possible remaining indflux
    idx_mat = false(size(ker2,1),num_ri-num_ri_now);
    for i=1:num_ri-num_ri_now
       idx_tmp = any(ker2(:,idx_list_rem(i)),2); 
       idx_mat(:,i) = idx_tmp;
    end

    idx_list = find(any(idx_mat,2));
    % X.rxn.rxn_names(idx_list);
    [idx_list_included,~] = ismember(idx_list,X.idx.idx_include);
    idx_list = idx_list(idx_list_included);

    b = nchoosek(idx_list,num_ri-num_ri_now);
    is_ok = false(1,size(b,1));
    for i=1:size(b,1)
       idx_now = idx_mat(b(i,:),:);
       is_ok(i) = all(any(idx_now,1));
    end
    b_ok = b(is_ok,:);

    % choose combination including the specified indflux
    num_ri_tmp = length(must_rxn)-num_ri_now;
    for j=1:num_ri_tmp
        idx_now = find(ismember(X.rxn.rxn_names,must_rxn{idx_must_rem(j)}));
        is_ok = false(1,size(b_ok,1));
        for i=1:size(b_ok,1)
           is_ok(i) = any(b_ok(i,:)==idx_now);
        end
        b_ok2 = b_ok(is_ok,:);
        b_ok = b_ok2;
    end

    b_ok_flux = arrayfun(@(x) X.rxn.rxn_names(x),b_ok,'UniformOutput',true);

    is_match = false(1,size(b_ok,1));
    for i=1:size(b_ok,1)
        is_match(i) = all(X.idx.idx_indflux(idx_list_rem)' == b_ok(i,:));
    end

    is_ok_indflux = any(is_match);

end