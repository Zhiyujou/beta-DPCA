function [opt_beta] = CV_beta(x, NumC, HM_delta, beta_seq, r, ru, q, k_fold)

% ==== Input: ====
% x: n x p data matrix
% NumC: no. machines
% HM_delta: regularization term
% beta_seq: vector of beta values
% r: true rank
% ru: upper bound of r
% q: sampling dimension 
% k_fold: K-fold cross-validation

% ==== Output: ====
% opt_beta: optimal value of beta




[n, p] = size(x);

if NumC<=k_fold
    k_fold = NumC;
end
NumC_tr = NumC - floor(NumC/k_fold);
NumC_va = floor(NumC/k_fold);


dist = zeros(length(beta_seq), k_fold);
qt = round(quantile(0:n, [0:k_fold]/k_fold));
for kk=1:k_fold
    ind = [1:n];
    ind = (ind>qt(kk)).*(ind<=qt(kk+1));
    x_tr = x(ind==0,:);   % training data
    x_va = x(ind==1,:);   % validation data

    n_tr = size(x_tr,1); 
    n_va = size(x_va,1); 
    nk_tr = floor(n_tr/NumC_tr)*ones(1, NumC_tr) + [ones(1, mod(n_tr, NumC_tr)), zeros(1, NumC_tr- mod(n_tr, NumC_tr))]; 
    nk_va = floor(n_va/NumC_va)*ones(1, NumC_va) + [ones(1, mod(n_va, NumC_va)), zeros(1, NumC_va- mod(n_va, NumC_va))]; 

    for ii=1:NumC_va
        count = sum(nk_va(1:ii-1));
        x_tmp_va = x_va(count+1:count+nk_va(ii),:);
        data_local_va = x_tmp_va; 
        [~, s_tmp_va, H_tmp_va] = svds(data_local_va, r, "largest");
        v_beta_va(:,:,ii) = H_tmp_va;     
    end

    [v_beta_tr, s_beta_tr] = beta_truncated(x_tr, NumC_tr, nk_tr, HM_delta, r, ru, q, beta_seq);    
  

    for ibeta = 1:length(beta_seq)
        for ii=1:NumC_va
            dist(ibeta, kk) = dist(ibeta, kk) + norm(( v_beta_tr(:, 1:r, ibeta)*v_beta_tr(:, 1:r, ibeta)'...
                                                     - v_beta_va(:,:,ii)*v_beta_va(:,:,ii)'), 'fro')^2/(NumC_va);
        end
    end
end
mean_dist = mean(dist, 2);

opt_beta_ind = find(mean_dist==min(mean_dist));
opt_beta = beta_seq(opt_beta_ind(1));    
end









