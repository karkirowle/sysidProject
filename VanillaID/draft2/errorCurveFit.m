% Housekeeping
clc;
clear;
close all;

% Creating fit

% fileList = {'C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_26_May_2018_14_43_37_50_200.mat' , ...,
%     'C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_28_May_2018_11_57_08_50_200.mat'};
%

SNR = [0.5,1,4,10,100];

figure;
hold on;



% Running index for cell
load('totalSNR');

measurementID = [];
trainSNR = [];
trainX = [];
trainY = [];
legendVector = [];
colors = {'r','g','b','c','m'};
figure;
for j = 1:length(SNR)
    % Conveniently then train-val-test for each state
    measurementID = [1:50, 1:50, 1:50];
    trainSNR = ones(1,150)*SNR(j);
    
    % A logarithmic predictor of measurement and bias based on SNR
    %     trainX = cat(2,trainX,measurementID);
    %     trainY = cat(1,trainY, [ squeeze(mean(totalMseMatrix(j,:,:,1),2)) ; ...,
    %         squeeze(mean(totalMseMatrix(j,:,:,2),2)); ...,
    %         squeeze(mean(totalMseMatrix(j,:,:,3),2))]);
    %
    trainX = measurementID;
    trainY = [ squeeze(mean(totalMseMatrix(j,:,:,1),2)) ; ...,
        squeeze(mean(totalMseMatrix(j,:,:,2),2)); ...,
        squeeze(mean(totalMseMatrix(j,:,:,3),2))];
    % Function model for fitting
    modelfun = @(b,x) 4 * log(x + 11) .^ b(1);
    %modelfun = @(b,x) exp(-b(1)*x + b(2));
    
    rng('default');
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'cauchy';
    opts.MaxIter = 1000;
    beta0 = 1;
    beta{j} = nlinfit(trainX',trainY,modelfun,beta0,opts);
    hold on;
    % Let's see how it looks
    h1 = plot(1:50, modelfun(beta{j}, trainX(:,1:50)'),colors{j},'LineStyle', ...,
        '--','LineWidth', 2);
    hold on;
    % Compare to original
    h2 = plot(1:50, squeeze(mean(totalMseMatrix(j,:,:,:),2)),colors{j}, 'LineWidth', 0.5);
    legendVector = [legendVector, h1(1), h2(1)];
    
end



legend(legendVector, {['Fit for SNR = ', num2str(SNR(1))],  ...,
    ['SNR = ', num2str(SNR(1))],  ...,
    ['Fit for SNR = ', num2str(SNR(2))],  ...,
    ['SNR = ', num2str(SNR(2))],  ...,
    ['Fit for SNR = ', num2str(SNR(3))],  ...,
    ['SNR = ', num2str(SNR(3))] ...,
    ['Fit for SNR = ', num2str(SNR(4))],  ...,
    ['SNR = ', num2str(SNR(4))] ...,
    ['Fit for SNR = ', num2str(SNR(5))],  ...,
    ['SNR = ', num2str(SNR(5))]
    }, 'FontSize', 12);

title('Logarithmic polynomial fit and RNMSE for different SNR [All states]' , ...,
    'FontSize', 16);

xlabel('# of measurements added (Fisher)', 'FontSize', 12);
ylabel('RNMSE', 'FontSize', 12);
xlim([1 50]);
ylim([0 1.5]);
run('figureFormatter');

% Sccater plot for polynomial values

figure;

plot(10*log10(SNR), [beta{1}, beta{2}, beta{3}, beta{4}, beta{5}],'red', ...,
    'Marker', 'o', 'LineWidth', 1.5);
xlabel('SNR (dB)', 'FontSize',12);
ylabel('value of fitted coefficient p (inverse rate of convergence)', 'FontSize', 12);
title('Tendency of faster convergence is reflected in coefficients', 'FontSize', 16);
run('figureFormatter');


% [B, fitInfo] = lasso(trainX', trainY, 'CV', 10);
% idxLambda1SE = fitInfo.Index1SE;
% coef = B(:,idxLambda1SE);
% coef0 = fitInfo.Intercept(idxLambda1SE);
%
%
%
% for j = 1:length(SNR)
%     id = (j-1)*150 + 1;
%     plot(1:50, modelfun(beta, trainX(:,id:id+49)'), ...,
%         'LineWidth', 1.5);
%     hold on;
%     plot(1:50, squeeze(mean(totalMseMatrix(j,:,:,3),2)), 'LineWidth', 1.5);
%     run('figureFormatter');
% end