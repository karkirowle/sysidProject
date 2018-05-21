% Housekeeping
clc;
clear;
close all;

%% Load checkpoint
addpath('checkpoints');
load('run_18-May-2018_40_1.mat')
rmpath('checkpoints');
load('colormapStore');

% TODO: Nonlinear assignments to color map for zero
% TODO: log scale or pruning the matrix?
%% Read matrix
for i=1:length(SNR)
    %     estimateMatrix = log10(abs(squeeze(estimate(i,1,1:40,:,1))));
    %     estimateMatrix(abs(estimateMatrix) == Inf) = 0;
    %     estimateMatrix = estimateMatrix;
    estimateMatrix = log10(abs(squeeze(estimate(i,1,1:40,:,1))));
    %     estimateMatrix(1:40,:) = estimateMatrix./norm(estimateMatrix(1:40,:));
    %     estimateMatrix = abs(5*log10(abs(estimateMatrix)));
    %     estimateMatrix(abs(estimateMatrix) == Inf) = 0;
    figure;
    
    imagesc(estimateMatrix);
    xlabel('dictionary function', 'FontSize', 16);
    ylabel('number of added measurements (Fisher)', 'FontSize', 16);
    title(['Evolution of identified weights while adding mo re measurements SNR=', num2str(SNR(i))], ...,
        'FontSize', 20);
    run('figureFormatter');
    grid off;
    colormap(mycmap);
    caxis([-20 10])
    colorbar;
    
    % Explicit color assignments
    colorMatrix = zeros(size(estimateMatrix,1),size(estimateMatrix,2),3);
    [x,y] = ind2sub(size(estimateMatrix),find(estimateMatrix == 0));
    
    for j=1:length(x)
        colorMatrix(x(j),y(j),1) = 1;
        % disp([num2str(x(j)), '_', num2str(y(k))]);
    end
    figure;
    surf(1:27,1:40,estimateMatrix,colorMatrix,'LineStyle', 'none');
    set(gca,'Ydir','reverse')
    axis image;
    view(2);
    [X, map] = rgb2ind(colorMatrix,1024);
    colormap(map);
    colorbar;
    
end
