% Unpruned weights

% Housekeeping
clc;
clear;
close all;

% Load the unpruned
load('zoomUnpruned.mat')

% How the unpruned weights move?
for i=1:81
   weights(i,:) = estimateUnprunedStruct{8,i}; 
end
figure;
plot(1:81, log(weights));
xlim([0 10]);

% How the costs move?
for i=1:81
   costs(i,:) = costStruct{8,i}; 
end
figure;
plot(1:81, costs);

