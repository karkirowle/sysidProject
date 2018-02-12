% RungeKutta simulation example

% Housekeeping
clc;
clear;
close all;

% Parameters
numberOfGenes = 3;
tspan = [0 100];
initialConditions = [10,10,10];
noise = 3;

assert(length(initialConditions) == numberOfGenes, ...
              'You havent provided the right amount initial conditions for the simulation!');
fun1 = @(x) x;
fun2 = @(x) 1./(1+x);
basisSum{1} = fun1;
basisSum{2} = fun2;
geneIdentifier(1) = 1;
geneIdentifier(2) = 2;
% Ez már kapásból nem jó, mert nem tudhatod hogy a function handle ugyanazt
% kapja e, tehát a summingot már csak a függvényen belül tudod megcsinálni

[t,x]=ode45(@(t,x) odefun(t,x,noise, basisSum, geneIdentifier, numberOfGenes),tspan,initialConditions,noise);

% Aim: sum the basis functions which could depend on each other 
function [f ] = odefun(t,x,noise, basisFunctions, geneIdentifier,numberOfGenes)
    % Preallocation required for ode45 helper function
    f = zeros(numberOfGenes,1);
    for j=1:numberOfGenes
        for i=1:2
            f(j) = f(j) + basisFunctions{i}(x(geneIdentifier(i)));
        end
        f(j) = noise*randn;
    end
    % TONOTE: Wiener considerations have not been made here
end
