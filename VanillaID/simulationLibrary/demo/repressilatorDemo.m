% Repressilator

% Housekeeping
clc;
clear;
close all;

node = 5;
sim = geneGraph(node);

k = 40;
d = -1;
sim = sim.repression(1, node, k, 4);
sim = sim.degradation(1,d);

for i=2:node
    sim = sim.repression(i, i-1, k, 4);
    sim = sim.degradation(i,d);
end



% Let's simulate the gene regulatory network
initialConditions = 1:1:node;
simulationInterval = 0:0.01:100;
processNoise = 10;

[derivativeSeries, ...,
    timeSeries] = ...,
    sim.runRungeKutta(initialConditions, processNoise, simulationInterval);

figure;
plot(simulationInterval, timeSeries, 'LineWidth', 1.5)
run('figureFormatter');
