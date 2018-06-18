

% Housekeeping
clc;
clear;
close all;

% Load stuff
%load('C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_26_May_2018_14_43_37_50_200.mat')
%load('C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_26_May_2018_23_05_01_200_100.mat')
%load('C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_28_May_2018_11_57_08_50_200.mat')

fileList = {'C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_26_May_2018_14_43_37_50_200.mat' , ...,
    'C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_28_May_2018_11_57_08_50_200.mat'};
figure;
hold on;

measurements = 1:50;
r = 200;
numRealisations = 200;
% Running index for cell
k = 1;

for a=1:1
    load('eachSNRgf2');
    % Diff eq selector
    for i=2:2
        % Signal to noise ratio
        for j=1:length(totalSNR)
            if (j==5)
                allStd = std(squeeze(totalMse(j,1:r,:,i)));
            else
                allStd = std(squeeze(totalMse(j,:,:,i)));
            end
            idx = measurements;
            idx(2*j:10:idx(end)) = [];
            allStd(idx) = 0;
            
            %stdIdx(1:10:50) = 1:10:50;
            
            errorbar(measurements, squeeze(mean(totalMse(j,:,:,i),2)),  ...,
                allStd, 'LineWidth', 1.5);
            
            legendString{k} = ['SNR = ', num2str(totalSNR(j))];
            k = k+1;
        end
        
        
        legend(legendString);
        xlabel('# of measurements added with Fisher algorithm', 'FontSize', 12);
        ylabel('RNMSE', 'FontSize', 12);
        title(['Means and standard deviations of RNMSE curves (', ...,
            num2str(numRealisations), ...,
            ' trials) [State ', num2str(i), ']'], 'FontSize', 16);
        run('figureFormatter');
        xlim([28 50]);
        ylim([0 2]);
    end
end