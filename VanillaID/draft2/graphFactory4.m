% Quick reproduction of all SNRs

% Housekeeping
clc;
clear;
close all,

% Load files
load('totalSNR');
load('oracleLinear');
measurements = 1:50;
numRealisations = 200;
SNR = [0.5, 1, 4, 10, 100];



% Diff eq selector
for i=1:3
    figure;
    hold on;
    % Signal to noise ratio
    for j=1:5
        
        allStd = std(squeeze(totalMseMatrix(j,:,:,i)));
        
        idx = measurements;
        idx(7*j:25:idx(end)) = [];
        allStd(idx) = NaN;
        
        %stdIdx(1:10:50) = 1:10:50;
        
        errorbar(measurements, squeeze(mean(totalMseMatrix(j,:,:,i),2)),  ...,
            allStd, 'LineWidth', 1.5);
        
        
        
    end
    
    % Oracle STD interval
    oracleStd = squeeze(std(RNMSE(:,:,i),1));
    idx = measurements;
    idx(3:15:end) = [];
    oracleStd(idx) = NaN;
    
    % Linear Oracle bounds
    errorbar(measurements, squeeze(mean(RNMSE(:,:,i),1)), ...,
        oracleStd, '--', 'LineWidth', 1.5);
    legend({'SNR = 0.1', 'SNR = 1', 'SNR = 10', 'SNR = 4', 'SNR = 100', ...,
        'Oracle linear estimator at SNR = 100'}, 'FontSize', 12);
    xlabel('# of measurements added with Fisher algorithm', 'FontSize', 12);
    ylabel('RNMSE', 'FontSize', 12);
    title(['Means and standard deviations of RNMSE curves (', ...,
        num2str(numRealisations), ...,
        ' trials) [State ', num2str(i), ']'], 'FontSize', 16);
    ylim([0 2]);
    xlim([1 50]);
    run('figureFormatter');
end
