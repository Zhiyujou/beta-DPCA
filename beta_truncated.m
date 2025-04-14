function [eigvec, eigenval]= beta_truncated(x, NumC, nk, HM_delta, r, ru, q, beta_seq)
% Method: beta-DPCA

% ==== Input: ====
% x: n x p data matrix
% NumC: no. machines
% nk: NumC-vector with subsample size for each machine
% HM_delta: regularization term
% r: true rank
% ru: upper bound of r (take the leading-ru eigenvectors)
% q: sampling dimension 
% beta_seq: vector of beta values

% ==== Output: ====
% eigvec: p x ru x length(beta_seq)
% eigenval: ru x ru x length(beta_seq)


[n, p] = size(x);
beta_avr = zeros(p, p, length(beta_seq));
for ii=1:NumC
    count = sum(nk(1:ii-1));
    x_tmp = x(count+1:count+nk(ii),:);
    data_local = x_tmp; 
    [~, s_tmp, H_tmp] = svds(data_local, q, "largest",'MaxIterations', 500);

    for ibeta = 1:length(beta_seq)
        if beta_seq(ibeta)==0 
            % beta → 0; vN (geometric mean, exp-log)
            beta_delta = 0; 
            beta_avr(:,:,ibeta) = beta_avr(:,:,ibeta)...
                + H_tmp*logm(s_tmp^2 + beta_delta*eye(length(diag(s_tmp))))*H_tmp'/NumC;  
        elseif beta_seq(ibeta) >0
            % positive beta
            beta_delta = 0;
            beta_avr(:,:,ibeta) = beta_avr(:,:,ibeta)...
                + (H_tmp*(s_tmp^2 + beta_delta*eye(length(diag(s_tmp))))^(beta_seq(ibeta))*H_tmp')/NumC;     
        else
            % negative beta
            beta_delta = HM_delta;
            beta_avr(:,:,ibeta) = beta_avr(:,:,ibeta)...
                + (H_tmp*(s_tmp)^2* H_tmp' + beta_delta*eye(p))^beta_seq(ibeta)/NumC;  
        end
    end
end


for ibeta = 1:length(beta_seq)   
    if beta_seq(ibeta)==0 
        % beta → 0; vN
        [vv_beta, ss_beta] = eigs(expm(beta_avr(:,:,ibeta)), ru, "largestreal");   
        v_beta(:,:,ibeta) = vv_beta;
        s_beta(:,:,ibeta) = ss_beta;
    else
        [vv_beta, ss_beta] = svds(beta_avr(:,:,ibeta)^(1/beta_seq(ibeta)), ru, "largest");   
        v_beta(:,:,ibeta) = vv_beta;
        s_beta(:,:,ibeta) = ss_beta;
    end
end
eigvec = v_beta;
eigenval = s_beta;

end



