% Experiment 1 - Analytical approximation

% Housekeeping
clc;
clear;
close all;

% --------------- Two-gene repressing toy example from ppt ----------------

sim = geneGraph(2);
sim = sim.degradation(1, -1);
sim = sim.degradation(2, -1);
sim = sim.repression(1,2,1,2);
sim = sim.repression(2,1,1,2);
[derivativeSeries, timeSeries] = sim.runSimulation([10,10], 100, 0, 1);

figure;
plot(timeSeries(:,1), 'LineWidth', 1.5)

x = zeros(100,1);
x(1) = 10;
x(2) = 0;

for t=3:100
    x(t) = 1 - 1./(2+2*x(t-2).^2+x(t-2).^4);
end

hold on;
plot(x, 'LineWidth', 1.5);
title('Euler simulation produces similar results to analytical derivation')
xlabel('timesteps')
ylabel('Concentration of x_{2}(t)')
box off;