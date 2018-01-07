
% ---------------- This includes test for simulationGraph2 ----------------
clc;
clear;
close all;
sim = simGraph(1);

% Auto-repressor network with 1 gene
% funHandle = @(x) x;
% funHandle2 = @(x) 1./(1+x);
% sim = sim.addBasisFunction(funHandle, 1,1, -0.6);
% sim = sim.addBasisFunction(funHandle2, 1,1, 0.5);
% disp(sim.basisFunctions);
% disp(length(sim.basisFunctions)); 
% results = sim.runSimulation(0.1, 100, 0.001, 1);
% 
% figure(1);
% plot(results);

% Two gene repressing each other
% 
sim = simGraph(2);
degradationHandle = @(x) x;
hillHandle = @(x) 1./((1+x).^4);

sim = sim.addBasisFunction(degradationHandle,1,1,-0.2);
sim = sim.addBasisFunction(hillHandle,1,2, 0.5);
sim = sim.addBasisFunction(degradationHandle,2,2, -0.2);
sim = sim.addBasisFunction(hillHandle,2,1, 0.5);

Y = sim.runSimulation([20, 10], 100, 0, 2);
figure;
plot(Y(:,1));
figure;
plot(Y(:,2));