% Housekeeping
clc;
clear;
close all;

%% Load checkpoint
addpath('checkpoints');
%load('run_22_May_2018_20_19_31_100_1')
%load('checkpoints\run_22_May_2018_17_56_59_50_100.mat')
load('C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_25_May_2018_02_29_14_50_100.mat')
rmpath('checkpoints');
load('colormapStore3');

%% Read matrix
for i=1:length(SNR)
    for j=1:3
        %     estimateMatrix = log10(abs(squeeze(estimate(i,1,1:40,:,1))));
        %     estimateMatrix(abs(estimateMatrix) == Inf) = 0;
        %     estimateMatrix = estimateMatrix;
        estimateMatrix = log10(abs(squeeze(mean(estimate(i,1,:,:,j),1))));
        %estimateMatrix = abs(squeeze(mean(estimate(i,1,:,:,1),1)));
        
        %     estimateMatrix(1:40,:) = estimateMatrix./norm(estimateMatrix(1:40,:));
        %     estimateMatrix = abs(5*log10(abs(estimateMatrix)));
        %     estimateMatrix(abs(estimateMatrix) == Inf) = 0;
        figure;
        
        imagesc(estimateMatrix);
        xlabel('dictionary function', 'FontSize', 16);
        ylabel('number of added measurements (Fisher)', 'FontSize', 16);
        title(['Evolution of identified weights with new measurements SNR=', num2str(SNR(i)), ' [State ', num2str(j) ']'], ...,
            'FontSize', 20);
        run('figureFormatter');
        grid off;
        colormap(mycmap);
        caxis([-20 10])
        colorbar;
    end
    % Explicit color assignments
    %     colorMatrix = zeros(size(estimateMatrix,1),size(estimateMatrix,2),3);
    %     [x,y] = ind2sub(size(estimateMatrix),find(estimateMatrix == 0));
    %
    %     for j=1:length(x)
    %         colorMatrix(x(j),y(j),1) = 1;
    %         % disp([num2str(x(j)), '_', num2str(y(k))]);
    %     end
    %     figure;
    %     surf(1:27,1:40,estimateMatrix,colorMatrix,'LineStyle', 'none');
    %     set(gca,'Ydir','reverse')
    %     axis image;
    %     view(2);
    %     [X, map] = rgb2ind(colorMatrix,1024);
    %     colormap(map);
    %     colorbar;
    %
end
