% Housekeeping
clc;
clear;
close all;

h = figure;
set(gcf,'color','w');

filename = 'sparseanimation3.gif';
%load('C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_25_May_2018_02_29_14_50_100.mat')
load('checkpoints/run_28_May_2018_11_57_08_50_200');
%Load stuff

diffeqid = 2;
snrid = 1;
lambdaError = NaN(1,1001);
%lambdaError(1:100:1001) = 3.2785;

% Textbox preformatting
dim = [.15 .6 .3 .3];
p = annotation('textbox',dim,'String','Sparsity count: ~' ...,
    ,'FitBoxToText','on');
p.LineStyle = 'none';

interval = 1:5:50;

for i=interval
    plot(1:1001,derivativeSeries(:,diffeqid));
    hold on;
    scatter(idx(1:i), corrDer(idx(1:i),diffeqid),40,'r','square','filled','LineWidth', ...,
        1.5);
    weights = squeeze(estimate(snrid,1,i,:,diffeqid));
    sparsity = length(find(weights ~= 0));
    p.String = ['Sparsity count: ', num2str(sparsity)];
    errorbar(1:1001,Phi*weights, lambdaError, ...,
        'LineWidth', 1.5);
    title('Change of fit with new measurements', 'FontSize', 16);
    xlabel('time (s)');
    ylabel('Derivative of concentration');
    hold off;
    
    drawnow;
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.5);
    end
end