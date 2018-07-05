% Housekeeping
clc;
clear;
close all;

% Load file
load('C:\Users\Lenovo\Documents\Bencur\imperial\fourth_year\projekt\sysidProject\VanillaID\checkpoints\run_29_May_2018_18_08_40_100_100.mat')


figure;
hold on;
for b=6:10
    tempMseMatrix = mseMatrixStruct{b};
    plot(1:100,squeeze(mean(tempMseMatrix(:,:,1),1)), 'LineWidth', 1.5);
end

xlabel('# of measurements added (Fisher)', 'FontSize', 12);
ylabel('RNMSE', 'FontSize', 12);
legend('Order 6', 'Order 7', 'Order 8', 'Order 9', 'Order 10');
title('No significant change in RNMSE when adding new functions', 'FontSize', 16);
run('figureFormatter');