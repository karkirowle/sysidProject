% Experiment 2 - Motif identification on two-gene net within new simulation
% environment

% Housekeeping
clc;
clear;
close all;

% --------------- Two-gene repressing toy example -------------------------

sim = geneGraph(2);
sim = sim.degradation(1, -0.2);
sim = sim.degradation(2, -0.2);
sim = sim.repression(1,2,0.5,4);
sim = sim.repression(2,1,0.5,4);

% --------------- Two-gene repressing toy example -------------------------
[derivativeSeries, timeSeries] = sim.runSimulation([20,10], 100, 0, 1);

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


[interpet, estimate, cost] = interpret.reconstruct(timeSeries(1:99,1), derivativeSeries(:,1));
motifSeries = interpret.motifCalculation(0:0.1:100, estimate,1);
figure;
plot(motifSeries);
xlabel('concentration of Gene 1')
ylabel('concentration derivative effect on Gene 1')