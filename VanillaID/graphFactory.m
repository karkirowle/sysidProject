% Graph factory for RMSE SNR curves

%Housekeeping
clc;
clear;
close all;

load('checkpoints\run_22_May_2018_17_56_59_50_100.mat')

measurements =  1:50;
for i=1:3
    figure;
    for j=1:5
        hold on;
        % This part is responsible for the selection of the error bars so
        % that they are evenly spaced out
        
        allStd = std(squeeze(mseMatrix(j,:,:,i)));
        idx = (1:50);
        idx(4*j:20:50) = [];
        allStd(idx) = 0;
        errorbar(measurements,mean(squeeze(mseMatrix(j,:,:,i))), ...,
            allStd, 'LineWidth', 1.5);
    end
    legend({'SNR=50', 'SNR=40', 'SNR=30', 'SNR=20', 'SNR=10'}, 'FontSize', 16);
    
    xlabel('# of measurements added with Fisher algorithm', 'FontSize', 12);
    ylabel('RNMSE', 'FontSize', 12);
    title('Means and standard deviations of RMSE curves (100 realisations)', 'FontSize', 16);
    run('figureFormatter');
    % print(['apr20/diffeq' num2str(i)],'-dpng','-r0');
    %     figure; plot(1:10,mean(mseMatrix(1:5,:,i)), 'LineWidth', 1.5);
    %     ylim([0 2]);
    %     xlabel('# of measurements added with Fisher algorithm');
    %     ylabel('RNMSE');
    % print(['apr20/diffeqmean' num2str(i)],'-dpng','-r0');
end

%figure; hist(squeeze(mseMatrix(1,:,1,1)));