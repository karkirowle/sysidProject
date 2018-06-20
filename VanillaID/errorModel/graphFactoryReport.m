% Graph factory for inverse MAP maximiser

% Housekeeping
clc;
clear;
close all;

% Load
load('MAPdata.mat');

% Parameters
measurements = 1:50;
diffeq = 1;

% Figure 1 - Diagram for error curves



% Preallocation
legendString = cell(1, length(totalSNR));

for j=1:3
    figure;
    plot(1:50,squeeze(mean(totalMse(:,:,:,j),2)), 'LineWidth', 1.5);
    xlabel('number of measurements added (inverse MAP covariance)', 'FontSize', 16);
    ylabel('RNMSE', 'FontSize', 16);
    title(['RNMSE calculated from optimally sampled measurements [X_{', ...,
        num2str(j), '}]'], ...,
        'FontSize', 20);
    xlim([1 50]);
    legendCell = cellstr(num2str(totalSNR', 'SNR = %g'));
    legend(legendCell, 'FontSize', 16);
    run('figureFormatter');
    figureFileName = ['errorModel/results/RMSE_Results_State_', num2str(j)];
    
    saveas(gcf,figureFileName,'epsc')
end
%Figure 1b - Comparison of each curve with Oracle Lower Bound
for j=1:3
    for i=1:length(totalSNR)
        figure;
        oracleRNMSE = oracleCalculator(Phi, totalSNR(i), groundTruth(:,j), ...,
            derivativeSeries(:,j), 50, squeeze(allCorrDer(i,:,:,j)));
        plot(1:50,squeeze(mean(totalMse(i,:,:,j),2)), 'LineWidth', 1.5);
        hold on;
        plot(1:50,oracleRNMSE, 'LineWidth', 1.5);
        xlabel('number of measurements added (inverse MAP covariance)', 'FontSize', 16);
        ylabel('RNMSE', 'FontSize', 16);
        ylim([0 2]);
        if (i == 1)
            ylim([0 10]);
        end
        title(['RNMSE curves while adding measurements [SNR = ', ...,
            num2str(totalSNR(i)), '] [X_{', num2str(j), '}]'], ...,
            'FontSize', 20);
        legend({'SBL', 'Oracle estimator'}, 'FontSize', 16);
        xlim([1 50]);
        figureFileName = ['errorModel/results/RMSE_Oracle_', strrep(num2str(totalSNR(i)),'.','_'), '_State_', num2str(j)];
        run('figureFormatter');
        saveas(gcf,figureFileName,'epsc')
    end
end



% Figure 2 - Diagram for ensemble dictionary visualisation

load('colorMapStore3');
for j=1:3    
    for i=1:length(totalSNR)
        estimateMatrix = log10(abs(squeeze(mean(allEstimate(i,:,:,:,j),2))));
        estimateMatrix = squeeze(estimateMatrix);
        % Threshold for ensemble pattern
        estimateMatrix(estimateMatrix < -1) = -Inf;
        
        figure;
        imagesc(estimateMatrix);
        xlabel('weight associated to the columns of the dictionary', 'FontSize', 16);
        ylabel('number of added measurements (inverse cov)', 'FontSize', 16);
        title(['Evolution of identified weights with measurements SNR=', ...,
            num2str(totalSNR(i)), ' [X_{', num2str(j) '}]'], ...,
            'FontSize', 20);
        
        run('figureFormatter');
        grid off;
        colormap(mycmap);
        caxis([-20 3])
        colorbar;
        figureFileName = ['errorModel/results/RMSE_Dictionary_SNR_', strrep(num2str(totalSNR(i)),'.','_'), '_State_', num2str(j)];
        saveas(gcf,figureFileName,'epsc')
        
    end
end

% Figure 3 - Diagram for informational quality ?

figure;
plot(1:50, squeeze(log10(mean(allFisher(:,:,1:50),2))), 'LineWidth', 1.5);
xlabel('number of measurements added (inverse MAP covariance)', 'FontSize', 16);
ylabel('Inverse MAP covariance', 'FontSize', 16);
title('Inverse covariance maximisation', ...,
    'FontSize', 20);
legendCell = cellstr(num2str(totalSNR', 'SNR = %g'));
legend(legendCell, 'FontSize', 16);
run('figureFormatter');

saveas(gcf,'errorModel/results/MAP_Covariance_SNR','epsc')