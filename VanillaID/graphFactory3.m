

% Housekeeping
clc;
clear;
close all;

% Load stuff
load('C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_25_May_2018_02_29_14_50_100.mat')


% Diff eq selector
for i=1:3
    figure;
    hold on;
    % Signal to noise ratio
    for j=1:4
        if (j==4)
            allStd = std(squeeze(mseMatrix(j,1:r,:,i)));
        else
            allStd = std(squeeze(mseMatrix(j,:,:,i)));
        end
        idx = (1:50);
        idx(1:10:50) = [];
        allStd(idx) = 0;
        
        stdIdx(1:10:50) = 1:10:50;
        if (i==4)
            errorbar(1:50, squeeze(mean(mseMatrix(j,1:r,:,i),2)),  ...,
                allStd, 'LineWidth', 1.5);
        else
            errorbar(1:50, squeeze(mean(mseMatrix(j,:,:,i),2)),  ...,
                allStd, 'LineWidth', 1.5);
        end
        legend({'SNR=50', 'SNR=40', 'SNR=30', 'SNR=20', 'SNR=10'}, 'FontSize', 16);
        xlabel('# of measurements added with Fisher algorithm', 'FontSize', 12);
        ylabel('RNMSE', 'FontSize', 12);
        title(['Means and standard deviations of RNMSE curves (100 trials) [State ', num2str(i), ']'], 'FontSize', 16);
        run('figureFormatter');
    end
end
