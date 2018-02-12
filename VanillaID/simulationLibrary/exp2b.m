% Experiment 2b - Motif identification on two-gene net within new simulation
% environment

% Housekeeping
clc;
clear;
close all;

% --------------- Two-gene repressing toy example -------------------------

% Parameters
initialConditions = [5,10];
timeStep = 1000;
precision = 0.1;
sim = geneGraph(2);
sim = sim.degradation(1, -0.8);
sim = sim.degradation(2, -0.8);
sim = sim.repression(1,2,10,2);
sim = sim.repression(2,1,10,2);

% --------------- Two-gene repressing toy example -------------------------
[derivativeSeries, timeSeries] = sim.runSimulation([40,10], timeStep, 0, precision);

% -------------- Interpretation graph test -------------------------------

interpret = interpretationGraph(1);
interpret = interpret.addBasisFunction(@(x) x);
interpret = interpret.addBasisFunction(@(x) 1./(1+x));
interpret = interpret.addBasisFunction(@(x) 1./(1+x).^2);
interpret = interpret.addBasisFunction(@(x) 1./(1+x).^3);
interpret = interpret.addBasisFunction(@(x) 1./(1+x).^4);
interpret = interpret.addBasisFunction(@(x) x./(1+x).^1);
interpret = interpret.addBasisFunction(@(x) x./(1+x).^2);
interpret = interpret.addBasisFunction(@(x) x./(1+x).^3);
interpret = interpret.addBasisFunction(@(x) x./(1+x).^4);


[interpet, estimate, cost] = interpret.reconstruct(timeSeries(1:(timeStep-1),1), derivativeSeries(:,1), 1);

% Motif plotted
motifSeries = interpret.motifCalculation(0:0.1:100, estimate,1);
figure(1);
plot(0:0.1:100,motifSeries, 'LineWidth', 1.5);
title('Motif')
xlabel('concentration of Gene 1')
ylabel('concentration derivative effect on Gene 1')
% Motif by measurements
motifSeries2 = interpret.motifCalculation(timeSeries(1:timeStep,1).', estimate,1);
figure(2);
plot(timeSeries(1:timeStep,1),motifSeries2, 'LineWidth', 1.5);
xlabel('concentration gene in measurement')
ylabel('respective derivative')
figure(3);
plot(1:timeStep,timeSeries(1:timeStep,1), 'LineWidth', 1.5);
xlabel('time')
ylabel('concentration')
title('Timeseries')
figure(4);
plot(1:timeStep-1,derivativeSeries(:,1), 'LineWidth',1.5);
hold on;
% plot(1:timeStep,motifSeries, 'LineWidth', 1.5);
% 
% title('Derivative');
