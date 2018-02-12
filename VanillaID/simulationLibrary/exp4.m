% Experiment 4 - Auto activation

% Housekeeping
clc;
clear;
close all;

% -- Auto activation toy example -- 

% Parameters
initialConditions = 10;
timeStep = 1000;
precision = 0.1;
sim = geneGraph(1);
sim = sim.degradation(1, -0.4);
sim = sim.activation(1,1,10,1);

% --------------- Two-gene repressing toy example -------------------------
[derivativeSeries, timeSeries] = sim.runSimulation(initialConditions, timeStep, 0, precision);

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


[interpet, estimate, cost] = interpret.reconstruct(timeSeries(1:(timeStep-1),1), derivativeSeries(:,1));

% Motif plotted
motifSeries = interpret.motifCalculation(0:0.1:100, estimate,1);
plot(timeSeries)