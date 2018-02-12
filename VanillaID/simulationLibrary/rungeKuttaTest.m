
%% ---------------- Runge Kutta testing ----------------

% Housekeeping
clc;
clear;
close all;

% Simple repressilator s
sim = simGraph(3);
degradationHandle = @(x) x;
hillHandle = @(x) 1./((1+x).^4);

sim = sim.addBasisFunction(degradationHandle,1,1,-0.1);
sim = sim.addBasisFunction(hillHandle,1,2, 0.7);
sim = sim.addBasisFunction(degradationHandle,2,2, -0.1);
sim = sim.addBasisFunction(hillHandle,2,3, 0.7);
sim = sim.addBasisFunction(degradationHandle,3,3, -0.1);
sim = sim.addBasisFunction(hillHandle,3,1, 0.7);

[derivativeSeries, timeSeries, timePoints] = sim.runRungeKutta([4, 2, 1], 0, [0 200]);
figure;
plot(timePoints, timeSeries);
figure;
plot(timePoints, derivativeSeries);