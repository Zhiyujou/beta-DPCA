%% Code for Simulation Studies in Section 4 %%

close all
clear all

rep = 100;
n = 250;   
p = 1000;
NumC = 5;  % no. machines
nk = floor(n/NumC)*ones(1, NumC) + [ones(1, mod(n, NumC)), zeros(1, NumC- mod(n, NumC))];  % subsample size for each machine

r = 5;  % true rank
q = r + 5;  % dimension of oversampling  (q = r+ r_os)
ru = r + 10;  %  upper bound of r (take the leading-ru eigenvectors)

HM_delta = 10^-5;  % regularization term
beta_seq = [-1,  0,  1];  % beta values                    
k_fold = 5;  % K-fold cross-validation

dof = 300;  % degrees of freedom for data distribution  

%%
c = p/n;
Lambda_critical = 1*(1+c^0.5);
Lambda_signal = Lambda_critical + p.^(1./[2:r+1])';
theta = 0.5*1;   % noise eigenvalues ~ U(theta, 2-theta)
Lambda = [Lambda_signal; unifrnd(theta, 2-theta, p-length(Lambda_signal),1)];

[Gamma, ~] = qr(normrnd(0,1,p,p));
Sigma = Gamma*diag(Lambda)*Gamma';
Sigma_half = Gamma*diag(Lambda.^(0.5))*Gamma';

%%
for irep = 1:rep
    % generate data    
    x = mvtrnd(eye(p), dof, n)*Sigma_half./(dof/(dof-2))^0.5;  
    
    
    %% DPCA Methods ( beta-DPCA,  Fan )
    % beta-DPCA
    for ibeta = 1:length(beta_seq)
        tic
        [v_beta_tmp, s_beta_tmp] = beta_truncated(x, NumC, nk, HM_delta, r, ru, q, beta_seq(ibeta));     
        time_beta(irep, ibeta) = toc;
    
        v_beta(:,:,ibeta) = v_beta_tmp;
        s_beta(:,:,ibeta) = s_beta_tmp;
    end
    
    % cv beta 
    tic
    [opt_beta] = CV_beta(x, NumC, HM_delta, beta_seq, r, ru, q, k_fold);
    [v_beta_opt, s_beta_opt] = beta_truncated(x, NumC, nk, HM_delta, r, ru, q, opt_beta);     
    time_beta_cv(irep) = toc;
    opt_beta_values(irep,:) = opt_beta;  

    
    % Fan
    tic
    [vv_f, ss_f] = Fan(x, NumC, nk, r, ru);
    time_fan(irep) = toc;

    %%
    for jj = r:ru
        for ibeta = 1:length(beta_seq)
            mean_similarity(irep,ibeta,jj) = mean(svds(v_beta(:,1:jj,ibeta)'*Gamma(:,1:r), r));  
        end
        mean_similarity(irep,length(beta_seq)+1,jj) = mean(svds(v_beta_opt(:,1:jj)'*Gamma(:,1:r), r));  
        mean_similarity(irep,length(beta_seq)+2,jj) = mean(svds(vv_f(:,1:jj)'*Gamma(:,1:r), r));  % fan (r)
    end
    irep
end

fprintf(' ======== freq beta =========\n')
tabulate(opt_beta_values(:,1))

auto_var_beta = arrayfun(@(x) ['beta = ', num2str(beta_seq(x))], 1:length(beta_seq), 'UniformOutput', false); 
disp(' * Average time (in sec.):')
disp(array2table([mean(time_fan), mean(time_beta,1), mean(time_beta_cv)], 'VariableNames', [{'Fan'}, auto_var_beta, {'beta (cv)'}]));

%%
figure;
hold on;
box on;
Mean_similarity = squeeze(mean(mean_similarity,1));

h1 = plot([r:ru], Mean_similarity(1,r:ru), '-.*b');  % beta-DPCA
h2 = plot([r:ru], Mean_similarity(2,r:ru), '-.ob'); 
h3 = plot([r:ru], Mean_similarity(3,r:ru), '-.xb');
h4 = plot([r:ru], Mean_similarity(4,r:ru), '-.dr','LineWidth',1);  % cv beta 
h5 = plot([r:ru], Mean_similarity(5,r:ru), ':ks','LineWidth',1);  % fan (r)  

legend([h5, h1, h2, h3, h4], 'Fan',...
    ['\beta = ',num2str(beta_seq(1))], ...
    ['\beta → ',num2str(beta_seq(2))], ...
    ['\beta = ',num2str(beta_seq(3))], ......
    '\beta_{cv}', ...
    'Location', 'SouthEast');
xlim([r-1, ru+1]);
ylim([0.0, 1.05]);
ylabel('\rho_k');
xlabel('k');

if dof==300
    title(['Gaussian dist.: (n, p, m) = (',num2str(n),', ',num2str(p),', ',num2str(NumC),')']);
else
    title(['t dist. with ', num2str(dof), ' dof: (n, p, m) = (',num2str(n),', ',num2str(p),', ',num2str(NumC),')']);
end


