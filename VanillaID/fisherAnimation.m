% Housekeeping
clc;
clear;
close all;

figure;
load('C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_25_May_2018_02_29_14_50_100.mat')

%Load stuff

diffeqid = 2;
snrid = 4;
lambdaError = NaN(1,1001);
lambdaError(1:100:1001) = lambda;

% Textbox preformatting 
dim = [.2 .8 .3 .3];
p = annotation('textbox',dim,'String',['Sparsity count: ~'] ...,
    ,'FitBoxToText','on');
p.LineStyle = 'none';

for i=1:size(estimate,3)
    plot(1:1001,derivativeSeries(:,diffeqid));
    hold on;
    scatter(idx(1:i), corrDer(idx(1:i),diffeqid));
    weights = squeeze(estimate(snrid,1,i,:,diffeqid));
    sparsity = length(find(weights ~= 0));
    p.String = ['Sparsity count: ', num2str(sparsity)];
    errorbar(1:1001,Phi*weights, lambdaError, ...,
        'LineWidth', 1.5);
    title('Change of fit with new measurements', 'FontSize', 16);
    xlabel('time (s)');
    ylabel('Derivative of concentration');
    hold off;
    pause(1);
end