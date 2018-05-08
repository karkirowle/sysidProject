% The vecrtor norms effect 2

% Housekeeping
clc;
clear;
close all;

load('zoomUnpruned');

for i=1:81
    % Question 1 - How the norm changes overall?
    estimateTemp2(i,:) = estimateUnprunedStruct{8,i};
    estimatenorm2(i) = norm(estimateTemp2(i,:));
end

figure;
plot(1:81,estimatenorm2);