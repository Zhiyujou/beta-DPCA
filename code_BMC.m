clear all
close all


load video_still_cell.mat  % 1-10 sec.
load video_data_cell.mat   % 11-20 sec.


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
frame_original = frame_original.*255;     
Y4beta = vertcat(Y{:}).*255;

Y_still = video_data_still_cell;
stillMean = zeros(480,640);
for i = 1:100
    tmp = Y_still{i}';
    stillMean = stillMean+tmp/100;
end
stillMean = stillMean.*255;    

NumC = size(Y,2);
nk = size(Y{1,1},1)*ones(1, NumC);   
p = size(Y{1,1},2);

%% plot original/background/moving part frames
%{/
figure;
tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'tight');
nexttile;
imshow(uint8(frame_original));
title(['Original frame'], 'FontSize', 14);

nexttile;
imshow(uint8(stillMean));
title(['Background'], 'FontSize', 14);

nexttile;
res = frame_original-stillMean;
res = (res-min(res(:)))./range(res(:));
res = res.*255;
imshow(uint8(res));
title(['Moving objects'], 'FontSize', 14);

obj_row = [300:380];
obj_col = [150:250];
nexttile;
imshow(uint8(res(obj_row,obj_col)));
title(['Zoomed-in view of the moving object'], 'FontSize', 14);
set(gcf, 'Position', [50, 250, 1700, 550]);  
%}


%% beta-DPCA
for ibeta = 1:length(beta_seq)
    tic
    [U_avg_tmp, s_beta_tmp] = beta_truncated(Y4beta, NumC, nk, HM_delta, basis_seq, basis_seq, q, beta_seq(ibeta));  
    time_beta(ibeta) = toc;
    U_beta(:,:,ibeta) = U_avg_tmp;
    reconstruct_global(:,:,ibeta) = U_avg_tmp*U_avg_tmp'*frame_original;
    residual_obj(:,:,ibeta) = (frame_original - reconstruct_global(:,:,ibeta));

    residual_tmp = residual_obj(:,:,ibeta);
    residual_tmp = (residual_tmp-min(residual_tmp(:)))./range(residual_tmp(:));
    residual_tmp = residual_tmp.*255;
    recon_residual_error(:, ibeta) = norm(res(obj_row,obj_col) - residual_tmp(obj_row,obj_col),'fro')/sqrt((range(obj_row)+1)*(range(obj_col)+1));
end

%% Fan
tic
[U_fan, ss_f] = Fan(Y4beta, NumC, nk, basis_seq, basis_seq);
time_fan = toc;
reconstruct_global(:,:,length(beta_seq)+1) = U_fan(:,1:basis_seq)*U_fan(:,1:basis_seq)'*frame_original;
residual_obj(:,:,length(beta_seq)+1) = (frame_original - reconstruct_global(:,:,length(beta_seq)+1));
  
residual_tmp = residual_obj(:,:,length(beta_seq)+1);
residual_tmp = (residual_tmp-min(residual_tmp(:)))./range(residual_tmp(:));
residual_tmp = residual_tmp.*255;
recon_residual_error(:,length(beta_seq)+1) = norm(res(obj_row,obj_col) - residual_tmp(obj_row,obj_col),'fro')/sqrt((range(obj_row)+1)*(range(obj_col)+1));


%%
figure;
t = tiledlayout(3, length(beta_seq)+1, 'Padding', 'compact', 'TileSpacing', 'tight');
nexttile(1);
tmp = reconstruct_global(:,:,length(beta_seq)+1);   
tmp = (tmp-min(tmp(:)))./range(tmp(:));
tmp = tmp.*255;
imshow(uint8(tmp));
title(['Fan'], 'FontSize', 16);

nexttile(1+(length(beta_seq)+1));
residual_map = abs(residual_obj(:,:,length(beta_seq)+1)); 
tmp = residual_map;
imshow(uint8(tmp));
title(['Fan'], 'FontSize', 16);

nexttile(1+2*(length(beta_seq)+1));
imshow(uint8(tmp(obj_row,obj_col)));
title(['Fan'], 'FontSize', 16);
%{/ 
dim1 = [0.056 0.61 0.3 0.1];    
annotation('textbox', dim1+[0 -0.61 0 0], ...
           'String', ['Reconstruction error : ' ], ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'left','FontSize', 14);

dim2 = [0.2 0.61 0.3 0.1];  
annotation('textbox', dim2+[0 -0.61 0 0], ...
           'String', [num2str(recon_residual_error(:,length(beta_seq)+1), '%.3f')], ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'left','FontSize', 14);
%}

for ibeta = 1:length(beta_seq)
    nexttile(1+ibeta);
    tmp = reconstruct_global(:,:,ibeta);  
    tmp = (tmp-min(tmp(:)))./range(tmp(:));
    tmp = tmp.*255;
    imshow(uint8(tmp));
    if beta_seq(ibeta)==0
        title(['\beta → ',num2str(beta_seq(ibeta)),''], 'FontSize', 16);
    else
        title(['\beta = ',num2str(beta_seq(ibeta)),''], 'FontSize', 16);
    end

    nexttile(1+(length(beta_seq)+1)+ibeta);
    residual_map = abs(residual_obj(:,:,ibeta));   
    tmp = residual_map;
    imshow(uint8(tmp));
    if beta_seq(ibeta)==0
        title(['\beta → ',num2str(beta_seq(ibeta)),''], 'FontSize', 16);
    else
        title(['\beta = ',num2str(beta_seq(ibeta)),''], 'FontSize', 16);
    end

    nexttile(1+2*(length(beta_seq)+1)+ibeta);
    imshow(uint8(tmp(obj_row,obj_col)));
    if beta_seq(ibeta)==0
        title(['\beta → ',num2str(beta_seq(ibeta)),''], 'FontSize', 16);
    else
        title(['\beta = ',num2str(beta_seq(ibeta)),''], 'FontSize', 16);
    end
    hold on;
   
    %{/
    dim3 = [0.133+0.233*(ibeta) 0.61 0.3 0.1];  
    annotation('textbox', dim3+[0 -0.61 0 0], ...
            'String', [num2str(recon_residual_error(:, ibeta), '%.3f')], ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'left','FontSize', 14);
    %}
end
set(gcf, 'Position', [50, 10, 1300, 990]);  

title_1 = ['(a) Projection to the leading rank-', num2str(basis_seq), ' eigenspace'];
text(-120, -250, title_1, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

title_2 = ['(b) Residual map'];
text(-130, -130, title_2, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

title_3 = ['(c) Zoomed-in region of the moving object'];
text(-120, -20, title_3, 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

