
clear all
close all

load olivettifaces    
rng(14685);

%% set parameters
NumC = 10;
basis_seq = 30;   % basis
q = basis_seq + 5; 
HM_delta = 10^(-5);
beta_seq = [-1, 0, 1];

%% data
faces = reshape(faces, [64^2,400]);
idx = 63;  % choose a frame for plotting reconstructions

%% noisy images
n_noise = 20;  % no. heterogeneous outiers  
noise_level = 50;  
mvt_dof = 5; 
imgStk = zeros(64^2, n_noise); 
for i = 1:n_noise
    tmp = randsample(0:255, 64^2, true);
    imgStk(:,i) = tmp(:);
end
faces_original = faces;
faces = [faces, imgStk];

%%
[p, n] = size(faces);
nk = floor(n/NumC)*ones(1, NumC) + [ones(1, mod(n, NumC)), zeros(1, NumC- mod(n, NumC))];  % subsample size for each machines

nn = size(faces, 2);
ind = randperm(nn);
faces_X = faces(:, ind);

%% 
Vhat = zeros(p, 1+length(beta_seq)); 
%% PCA
[U_pca, ss_pca] = svds(cov(faces_X'), basis_seq, "largest");
Vhat(:, 1) = U_pca(:, 1:basis_seq)*(U_pca(:, 1:basis_seq)'*faces(:, idx));
recon_error(:, 1) = norm(faces(:, idx) - Vhat(:, 1),'fro')/sqrt(64^2);

%% beta-method
for ibeta = 1:length(beta_seq)
    tic
    [U_avg_tmp, s_beta_tmp] = beta_truncated(faces_X', NumC, nk, HM_delta, basis_seq, basis_seq, q, beta_seq(ibeta));    
    time_beta(ibeta) = toc;
    U_beta(:,:,ibeta) = U_avg_tmp;
    Vhat(:, ibeta+1) = U_avg_tmp(:, 1:basis_seq)*(U_avg_tmp(:, 1:basis_seq)'*faces(:, idx));
    recon_error(:, ibeta+1) = norm(faces(:, idx) - Vhat(:, ibeta+1),'fro')/sqrt(64^2);
end


%%
figure;
box on;
hold on;
h = 1; 
w = size(Vhat, 2)+1;
subplot(h, w, 1);

tmp = reshape(faces(:,idx), [64,64]);
imshow(uint8(tmp));
title(['True'],'FontSize', 12);
for i = 1:size(Vhat, 2)
    subplot(h, w, i+1);
    tmp = squeeze(Vhat(:, i));
    tmp = reshape(tmp, [64,64]);
    imshow(uint8(tmp));

    if i == 1
        title(['PCA'],'FontSize', 12);
    elseif i>1 
        if beta_seq(i-1) == 0 
            title(['\beta → ', num2str(beta_seq(i-1))], 'FontSize', 12);
        else
            title(['\beta = ', num2str(beta_seq(i-1))], 'FontSize', 12);
        end
    end
    dim1 = [0.13 0.1 0.3 0.1];   
    annotation('textbox', dim1, ...
            'String', ['Reconstruction error : ' ], ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'left', 'FontSize', 12);
    dim2 = [0.325+0.165*(i-1) 0.1 0.3 0.1];  
    annotation('textbox', dim2, ...
            'String', [num2str(recon_error(:, i), '%.3f')], ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'left', 'FontSize', 12);
end
stitle = sgtitle(['Using ', num2str(basis_seq), ' PCs'], 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
set(gcf, 'Position', [50, 200, 1300, 300]);   

