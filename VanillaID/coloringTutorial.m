% Housekeeping
clc;
clear;
close all;

%% Load checkpoint
addpath('checkpoints');
load('run_18-May-2018_40_1.mat')
rmpath('checkpoints');

%% Cucc
estimateMatrix = abs(squeeze(estimate(1,1,1:40,:,1)));

cMap = parula(256);

% Maximum value of data
dataMax = max(max(estimateMatrix));

% Minimum value of data
dataMin = min(min(estimateMatrix));

centerPoint = 10^(-10);
centerPoint2 = 1;


scalingIntensity = 4;

% Interval of X
x = 1:length(cMap);

% This

x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
 x = scalingIntensity * x/max(abs(x));

figure; plot(x);
% x = scalingIntensity * x/max(abs(x));
% 
% x = sign(x).* exp(abs(x));

x = x - min(x); x = x*511/max(x)+1;
figure; plot(x);

% Linearly interpolate the values on x to the colormap
newMap = interp1(x, cMap, 1:512);



figure; imagesc(estimateMatrix);

figure; imagesc(estimateMatrix); colormap(newMap);