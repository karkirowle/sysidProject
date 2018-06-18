% Housekeeping
clc;
clear;
close all;

% Load the merged experiment data
load('eachSNR');


legendVector = [];
legendCell = cell(1,length(totalSNR));
beta = zeros(length(totalSNR),2);
% Color vector containg the colors of the lines
colors = {'r','g','b','c','m','y','k','r','g','b','c'};
subselect = 1;
figure;
for j = 1:length(totalSNR)
    % Conveniently then train-val-test for each state
    
    trainX = repmat(subselect:50,1,3);
    
    trainY = [ squeeze(mean(totalMse(j,:,subselect:end,1),2)); ...,
        squeeze(mean(totalMse(j,:,subselect:end,2),2)); ...,
        squeeze(mean(totalMse(j,:,subselect:end,3),2))];
    
%         trainY = [ reshape(squeeze(totalMse(j,:,:,1)), [1 50*200]),  ...,
%             reshape(squeeze(totalMse(j,:,:,2)), [1 50*200]) ...,
%             reshape(squeeze(totalMse(j,:,:,3)), [1 50*200])];
    
    
    % Function model for fitting
    modelfun = @(b,x) b(2) * log(x.^2.3 + 11) .^ b(1);
    %modelfun = @(b,x) b(1).*log(x + b(3)).^(b(2));
    %modelfun = @(b,x) b(1).*squeeze(mean(totalMse(1,:,x,1),2))' + b(2);
    %modelfun = @(b,x) b(2)*(1 - b(1)).^(2.*x) + b(2) ;
%     modelfun = @(b,x) b(1).*(1 - b(2)).^(2.*x) + b(3);

    rng('default');
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'cauchy';
    opts.MaxIter = 1000;
    beta0 = [-1,1];
    
    [beta(j,:),~,J,~,MSE(j),~] = nlinfit(trainX',trainY,modelfun,beta0,opts);
    hold on;
    
    % Let's see how it looks
    h1 = plot(1:50, modelfun(beta(j,:), trainX(:,subselect:50)'),colors{j},'LineStyle', ...,
        '--','LineWidth', 2);
    hold on;
    % Compare to original states
    h2 = plot(1:50, squeeze(mean(totalMse(j,:,:,:),2)),colors{j}, ...,
        'LineWidth', 0.5);
    
    % Save legend object references
    legendVector = [legendVector, h1(1), h2(1)];
    
    % Save legend strings
    legendCell{2*j-1} = ['Fit for SNR = ', num2str(totalSNR(j))];
    legendCell{2*j} = ['SNR = ', num2str(totalSNR(j))];
    
end


legend(legendVector, legendCell, 'FontSize', 12);
title('Logarithmic polynomial fit and RNMSE for different SNR [All states]' , ...,
    'FontSize', 16);
xlabel('# of measurements added (Fisher)', 'FontSize', 12);
ylabel('RNMSE', 'FontSize', 12);
xlim([1 50]);
run('figureFormatter');

% Sccater plot for polynomial values

figure;
plot(10*log10(totalSNR), real(beta(:,1)) ,'red', ...,
    'Marker', 'o', 'LineWidth', 1.5);
xlabel('SNR (dB)', 'FontSize',12);
ylabel('value of fitted coefficient p (inverse rate of convergence)', 'FontSize', 12);
title('Tendency of faster convergence is reflected in coefficients', 'FontSize', 16);
run('figureFormatter');

figure;
plot(10*log10(totalSNR), real(beta(:,2)) ,'red', ...,
    'Marker', 'o', 'LineWidth', 1.5);
xlabel('SNR (dB)', 'FontSize',12);
ylabel('value of fitted coefficient p (inverse rate of convergence)', 'FontSize', 12);
title('Tendency of faster convergence is reflected in coefficients', 'FontSize', 16);
run('figureFormatter');

figure;
plot(10*log10(totalSNR),MSE); 
