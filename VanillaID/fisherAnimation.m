% Bayesian Identification of Gene Regulatory Networks
% Bence Halpern 2018

% Housekeeping
clc;
clear;
close all;


h = figure;
% Sets the background of the figure white
set(gcf,'color','w');

% Loading parameters
filename = 'sparseanimation3.gif';
load('checkpoints/run_28_May_2018_11_57_08_50_200');
load('utilities/colorMapStore3.mat');

% Parameters which tell the program what to display

% This selects the gene (1-3)
diffeqid = 2;

% This selects the SNR (1-length(SNR))
snrid = 1;

% TODO: Why NAN array here?
lambdaError = NaN(1,1001);
%lambdaError(1:100:1001) = 3.2785;

% Place the textbox on the upper left corner of the figure 

dim = [.15 .6 .3 .3];
p = annotation('textbox',dim,'String','Sparsity count: ~' ...,
    ,'FitBoxToText','on');
p.LineStyle = 'none';

% Measurement interval (TODO: Not necessarily best way to do this)
interval = 1:1:50;

% TODO: again, this is hard-coded, not a nice way
weightMatrix = ones(50,27) * -Inf;
% We are creating two subplots and update it, and then append it to
% the total gif

for i=interval
    
    % Left plot - How the fit changes
    subplot(1,2,1);
    
    % Show the ground truth for a given differential equation
    plot(1:1001,derivativeSeries(:,diffeqid));
    hold on;
    
    % Show what are the sampled noisy values
    scatter(idx(1:i), corrDer(idx(1:i),diffeqid),40,'r','square','filled','LineWidth', ...,
        1.5);
        
    % Fetch the weights
    weights = squeeze(estimate(snrid,1,i,:,diffeqid));
    
    % Update the textbox with the number of non-zero elements
    sparsity = length(find(weights ~= 0));
    p.String = ['Sparsity count: ', num2str(sparsity)];
    
    % And finally show the results of identification algorithm
    errorbar(1:1001,Phi*weights, lambdaError, ...,
        'LineWidth', 1.5);
    title('Change of fit with new measurements', 'FontSize', 16);
    xlabel('time (s)');
    ylabel('Derivative of concentration');
    hold off;
    
    % Right plot - How the weights evolve
    subplot(1,2,2);
    
    % Using log to shrink the magnitudes, and we repeatedly fetch
    % from 1:i TODO: could be simplfied to (i from 1:i)
    weightMatrix(1:i,:) = log10(abs(squeeze(estimate(snrid,1,1:i,:,diffeqid))));
    imagesc(weightMatrix);
    colormap(mycmap);
    title('Evolution of weights in dictionary', 'FontSize', 16);
    xlabel('dictionary function id');
    ylabel('# of measurements added (Fisher)');
    set(gcf,'Position',[50 50 1500 1000]);
    drawnow;
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if i == 1
        % On first iteration just create a lWooping gif
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        % On other iterations apped the frame to this gif
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.5);
    end
end