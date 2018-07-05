% Bayesian Identification of Gene Regulatory Networks
% Bence Halpern 2018

% Housekeeping
clc;
clear;
close all;

%% Load checkpoint
addpath('checkpoints');

% Write the desired simulation name in load
load('run_28_May_2018_11_57_08_50_200');
rmpath('checkpoints');
load('colormapStore3');

%% Read matrix
for i=1:length(SNR)
    for j=1:3
        
        % Taking the log of the mean of the weights to compress
        % changes in scale
        estimateMatrix = log10(abs(squeeze(mean(estimate(1,:,:,:,j),2))));
        estimateMatrix = squeeze(estimateMatrix);
        
        % Threshold for ensemble pattern
        estimateMatrix(estimateMatrix < 0) = -Inf;
       
        % Code for generating the actual plot
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
end
