function [eigvec, eigenval]= Fan(x, NumC, nk, r, ru)
% Method: DPCA by Fan et al. (2019)

% ==== Input: ====
% x: n x p data matrix
% NumC: no. machines
% nk: NumC-vector with subsample size for each machine
% r: true rank
% ru: upper bound of r (take the leading-ru eigenvectors)

% ==== Output: ====
% eigvec: p x ru 
% eigenval: ru x ru 


[n, p] = size(x);
f_avr = zeros(p,p);
for ii=1:NumC
    count = sum(nk(1:ii-1));
    x_tmp = x(count+1:count+nk(ii),:);
    data_local = x_tmp;
    [~, s_tmp, H_tmp] = svds(data_local, r, "largest");
    f_avr = f_avr + H_tmp*H_tmp'/NumC;  

end
[vv_f, ss_f] = svds(f_avr, ru);  

eigvec = vv_f;
eigenval = ss_f;

end

