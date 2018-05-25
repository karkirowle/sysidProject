% Fisher information graph factory

% Graph factory for RMSE SNR curves

%Housekeeping
clc;
clear;
close all;

load('checkpoints\run_22_May_2018_17_56_59_50_100.mat')

measurements =  1:50;
figure; 
for i=1:5
    hold on;
    plot(1:50, squeeze(log10(fisherDetMatrix(i,1,1:50))), 'LineWidth', 1.5);
    legend({'SNR=50', 'SNR=40', 'SNR=30', 'SNR=20', 'SNR=10'}, 'FontSize', 16);
    
    xlabel('# of measurements added with Fisher algorithm', 'FontSize', 12);
    ylabel('base 10 logarithm of Fisher Information', 'FontSize', 12);
    title('Fisher information monotone increases for each measurement point', 'FontSize', 16);
    run('figureFormatter');
    % print(['apr20/diffeq' num2str(i)],'-dpng','-r0');
    %     figure; plot(1:10,mean(mseMatrix(1:5,:,i)), 'LineWidth', 1.5);
    %     ylim([0 2]);
    %     xlabel('# of measurements added with Fisher algorithm');
    %     ylabel('RNMSE');
    % print(['apr20/diffeqmean' num2str(i)],'-dpng','-r0');
end

%figure; hist(squeeze(mseMatrix(1,:,1,1)));