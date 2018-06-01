% Oracle linear regression bounds

% Housekeeping
clc;
clear;
close all,

% Load data here
fileList = {'C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_26_May_2018_14_43_37_50_200.mat' , ...,
    'C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_28_May_2018_11_57_08_50_200.mat'};
load(fileList{2});

% Best estimator if subset is known
RNMSE = zeros(length(SNR),numRealisations, length(measurements), 3);
RNMSE2 = zeros(length(SNR), numRealisations, length(measurements), 3);
previous = [3, 1, 2];
for k = 1:length(SNR)
    for j = 1:3
        for r = 1:numRealisations
            for i=1:length(measurements)
                currentIdx = idx(1:i);
                corrDerAct = squeeze(corrDerMatrix(k,r,currentIdx,j));
                timeSeriesAct = timeSeries(currentIdx,j);
                timeSeriesAct2 = timeSeries(currentIdx,previous(j));
                % Degradation function (linear)
                degFun = (@(x) x);
                repFun = @(x) 1./(1+x.^4);
                
                Dic = [degFun(timeSeriesAct), repFun(timeSeriesAct2)];
                
                
                
                % Perform regression here
                estimateLSE = (Dic' * Dic) \ Dic' * corrDerAct;
                
                % Perform normal regression here
                estimateLSE2 = (Phi(1:i,:)' * Phi(1:i,:)) \ Phi(1:i,:)' * corrDerAct;
                
                % RMSE metric
                estimateLSEVector = zeros(27,1);
                estimateLSEVector(1:2) = estimateLSE;
                
                
                groundTruth = zeros(27,1);
                groundTruth(1) = -1;
                groundTruth(2) = 40;
                
                groundTruth2 = zeros(27,1);
                groundTruth2(j) = -1;
                groundTruth2(3*4 + j) = 40;
                
                RNMSE(k,r,i,j) = norm(estimateLSEVector-groundTruth,2)/norm(groundTruth,2);
                RNMSE2(k,r,i,j) = norm(estimateLSE2-groundTruth2,2)/norm(groundTruth2,2);
                
            end
        end
        figure;
        %plot(measurements, squeeze(mean(RNMSE2(:,:,j),1)), 'LineWidth', 1.5);
        hold on;
        plot(measurements, squeeze(mean(RNMSE(:,:,j),1)), 'LineWidth', 1.5);
        xlabel('# of measurements added');
        ylabel('RNMSE');
        run('figureFormatter');
    end
end
% Best estimator if subset is unknown



clearvars -except RNMSE RNMSE2
save oracleLinear

