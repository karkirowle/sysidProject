% Zooming into the weights

% Housekeeping
clc;
clear;
close all;

% Load data
load('zoomWeights2')
clear title;


weightVal1 = zeros(1,81);
weightVal2 = zeros(1,81);

for i=1:81
    estimateM = estimateStruct{i};
    weightVal1(i) = estimateM(4,4);
    weightVal2(i) = estimateM(8,8);
end

figure;
scatter(weightVal1, weightVal2);

xIndices = 1:10;
xIndices = [xIndices, 11:8:81];

x = weightVal1(xIndices(1:end-1));
y = weightVal2(xIndices(1:end-1));

scale = 1.5*(max(y)-min(y))/(max(x)-min(x));
x = scale*x;
u = scale*diff(weightVal1(xIndices(1:end)));
v = diff(weightVal2(xIndices(1:end)));


figure;
q = quiver(x,y,u,v, 'LineWidth', 1.5);
xlabel('$\hat{w_1}$ (697.5908 times upscaled) ','Interpreter','latex')
ylabel('$\hat{w_2}$','Interpreter','latex');
title('Convergence of weights at the cliff');
hold on;
scatter(scale*groundTruth(4,4),groundTruth(8,8),'filled');
legend('optimisation path', 'ground truth','Location','Northwest');

% Look at the weights

for i=1:81
    estimateM = estimateStruct{i};
    weightVals(1:72,i) = estimateM(:,8);
end

% 3,8, 31,35, 39
% 12 iterations

% figure;
% scatter(weightVals(3,1:12), weightVals(8,1:12));
% hold on;
% scatter(weightVals(31,1:12), weightVals(35,1:12));
% hold on;
% scatter(weightVals(39,1:12), weightVals(40,1:12));

scale = 12.5;
w3 = weightVals(3,1:11);
w8 = weightVals(8,1:11)*scale;
w31 = weightVals(31,1:11);
w35 = weightVals(35,1:11)*scale;
w39 = weightVals(39,1:11);
w40 = weightVals(40,1:11)*scale;


% Diff
dw3 = diff(weightVals(3,1:12));
dw8 = diff(weightVals(8,1:12))*scale;
dw31 = diff(weightVals(31,1:12));
dw35 = diff(weightVals(35,1:12))*scale;
dw39 = diff(weightVals(39,1:12));
dw40 = diff(weightVals(40,1:12))*scale;

figure;
quiver(w3,w8,dw3,dw8, 'LineWidth', 1.5);
hold on,
quiver(w31,w35,dw31,dw35, 'LineWidth', 1.5);
hold on,
quiver(w39,w40,dw39,dw40, 'LineWidth', 1.5);
hold on;
scatter(groundTruth(3,8),scale*groundTruth(8,8),'filled');
hold on;
scatter(groundTruth(31,8),scale*groundTruth(35,8),'filled');
hold on;
scatter(groundTruth(39,8),scale*groundTruth(40,8),'filled');
xlabel('$\hat{w_1}$ (697.5908 times upscaled) ','Interpreter','latex')
ylabel('$\hat{w_2}$','Interpreter','latex');
