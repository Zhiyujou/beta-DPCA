clear all
close all


load video_data_cell.mat

%% set parameters
basis_seq = 35;  % no. components
q = basis_seq + 25; 
HM_delta = 10^(-5);
beta_seq = [-1, 0, 1];

%% data
start_frame = 176;
Y = video_data_cell(start_frame:250);
idx = 40;  % choose a frame for plotting reconstructions
frame_original = Y{1,idx}';
Y4beta = vertcat(Y{:});

NumC = size(Y,2);
nk = size(Y{1,1},1)*ones(1, NumC);   % subsample size for each machine
p = size(Y{1,1},2);


%%
% beta-DPCA
for ibeta = 1:length(beta_seq)
    tic
    [U_avg_tmp, s_beta_tmp] = beta_truncated(Y4beta, NumC, nk, HM_delta, basis_seq, basis_seq, q, beta_seq(ibeta));  
    time_beta(ibeta) = toc
    U_beta(:,:,ibeta) = U_avg_tmp;
    reconstruct_global(:,:,ibeta) = U_avg_tmp*U_avg_tmp'*frame_original;
end

% Fan
tic
[U_fan, ss_f] = Fan(Y4beta, NumC, nk, basis_seq, basis_seq);
time_fan = toc
reconstruct_global(:,:,length(beta_seq)+1) = U_fan(:,1:basis_seq)*U_fan(:,1:basis_seq)'*frame_original;



%%
figure;
t = tiledlayout(2,4, 'Padding', 'compact', 'TileSpacing', 'tight');
nexttile(1)
tmp = reconstruct_global(:,:,4);   % background
imshow(tmp)
title(['Fan'], 'FontSize', 14)

nexttile(5)
residual_map = abs(double(frame_original) - double(tmp)); % Residual Map
imshow(residual_map)
title(['Fan'], 'FontSize', 14)


for jj=1:length(beta_seq)
    nexttile(1+jj);
    tmp = reconstruct_global(:,:,jj);  % background
    imshow(tmp)
    if beta_seq(jj)==0
        title(['\beta → ',num2str(beta_seq(jj)),''], 'FontSize', 14)
    else
        title(['\beta = ',num2str(beta_seq(jj)),''], 'FontSize', 14)
    end

    nexttile(5+jj);
    residual_map = abs(double(frame_original) - double(tmp));  % Residual Map
    imshow(residual_map)
    if beta_seq(jj)==0
        title(['\beta → ',num2str(beta_seq(jj)),''], 'FontSize', 14)
    else
        title(['\beta = ',num2str(beta_seq(jj)),''], 'FontSize', 14)
    end
end
title_1 = ['(a) Projection to the leading rank-', num2str(basis_seq), ' eigenspace'];
text(-700, -770, title_1, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

title_2 = ['(b) Residual map'];
text(-700, -100, title_2, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');








