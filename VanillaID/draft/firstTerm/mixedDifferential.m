% Consider a differential equation with activation and represison

% Housekeeping
clc;
clear;
close all;

% Diff simulation
interval = 0:0.001:10;

w1 = rand*100;
w2 = rand*100;
w3 = rand*100;
w4 = rand*100;
repression1 = @(x) 1./(1+x.^1);
activation1 = @(x) x.^1./(1+x.^1);
repression2 = @(x) 1./(1+x.^2);
activation2 = @(x) x.^2./(1+x.^2);
dxdt = w1.*repression1(interval) + w2.*activation1(interval) + w3.*repression2(interval) ...,
    +  w4.*activation2(interval);

figure(1);
plot(interval, dxdt);