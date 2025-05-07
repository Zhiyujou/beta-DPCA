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
obj_row = [300:380];
obj_col = [150:250];
% figure;
% imshow(frame_original(obj_row,obj_col));


figure;
t = tiledlayout(3,4, 'Padding', 'compact', 'TileSpacing', 'tight');
nexttile(1);
tmp = reconstruct_global(:,:,4);   % background
imshow(tmp);
title(['Fan'], 'FontSize', 16);

nexttile(5);
residual_map = abs(double(frame_original) - double(tmp)); % Residual Map
imshow(residual_map);
title(['Fan'], 'FontSize', 16);

nexttile(9);
imshow(residual_map(obj_row,obj_col));
title(['Fan'], 'FontSize', 16);


for jj = 1:length(beta_seq)
    nexttile(1+jj);
    tmp = reconstruct_global(:,:,jj);  % background
    imshow(tmp);
    if beta_seq(jj)==0
        title(['\beta → ',num2str(beta_seq(jj)),''], 'FontSize', 16);
    else
        title(['\beta = ',num2str(beta_seq(jj)),''], 'FontSize', 16);
    end

    nexttile(5+jj);
    residual_map = abs(double(frame_original) - double(tmp));  % Residual Map
    imshow(residual_map);
    if beta_seq(jj)==0
        title(['\beta → ',num2str(beta_seq(jj)),''], 'FontSize', 16);
    else
        title(['\beta = ',num2str(beta_seq(jj)),''], 'FontSize', 16);
    end
    nexttile(9+jj);
    imshow(residual_map(obj_row,obj_col))
    if beta_seq(jj)==0
        title(['\beta → ',num2str(beta_seq(jj)),''], 'FontSize', 16);
    else
        title(['\beta = ',num2str(beta_seq(jj)),''], 'FontSize', 16);
    end
    hold on;
   
end
set(gcf, 'Position', [50, 10, 1300, 990]);  

title_1 = ['(a) Projection to the leading rank-', num2str(basis_seq), ' eigenspace'];
text(-120, -250, title_1, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

title_2 = ['(b) Residual map'];
text(-130, -130, title_2, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

title_3 = ['(c) Zoomed-in region of the moving object'];
text(-120, -20, title_3, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');









