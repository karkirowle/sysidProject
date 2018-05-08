% The vecrtor norms effect

% Housekeeping
clc;
clear;
close all;

load('zoomWeights2');

for i=1:81
estimateTemp = estimateStruct{i};
estimatenorm(i) = norm(estimateTemp(:,8),2);
estimates(i,1) = estimateTemp(3,8)^2/estimatenorm(i)^2;
estimates(i,2) = estimateTemp(8,8)^2/estimatenorm(i)^2;
estimates(i,3) = estimateTemp(31,8)^2/estimatenorm(i)^2;
estimates(i,4) = estimateTemp(35,8)^2/estimatenorm(i)^2;
estimates(i,5) = estimateTemp(39,8)^2/estimatenorm(i)^2;
end

figure;
plot(1:81, estimatenorm, 'LineWidth', 1.5);
% hold on;
% plot(1:81, estimates(:,[4]), 'LineWidth', 1.5);
% hold on;
% plot(1:81, repmat(10e-06,[1,81]), 'LineWidth', 1.5);
% xlim([0 14]);