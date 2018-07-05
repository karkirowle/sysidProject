% Bayesian Identification of Gene Regulatory Networks

% Bence Halpern 2018

% ------------------------- DESCRIPTION -----------------------------------
% Fisher information based simulation algorithm of a repressilator network
% This file is for the parameter sweep of the noise
% Repressilator parameters are kept constant

% Housekeeping
clc;
clear;
close all;

% ------------------- PARAMETERS: CHANGE WISELY! --------------------------
measurements = 1:50;
SNR = [10,1,1/2];
numRealisations = 100;

% Server parameters
clusterNumber = 16;
parpool('local',clusterNumber)

% -------------------------------------------------------------------------

% Fetching date for filename
filenameDate = datetime('now','Format','default');
dateChar = char(filenameDate);
dateChar(dateChar == ' ') = '_';
dateChar(dateChar == '-') = '_';
dateChar(dateChar == ':') = '_';


% Setting the seed to ensure reproducibility of experimental results
rng('default');

% Creating three gene repressilator topology
nodes = 3;
sim = geneGraph(nodes);

% Don't change the one argument here
interpret = interpretationGraph(1);

% Adding basis functions to the repressilator model

% Degradation function (linear)
interpret = interpret.addBasisFunction(@(x) x);

% Hill functions added up to order 4
for i=1:4
    interpret = interpret.addBasisFunction(@(x) 1./(1+x.^i));
end

for i=1:4
    interpret = interpret.addBasisFunction(@(x) (x.^i)./(1+x.^i));
end

% Each gene has it's own basis function -> total number of functions
numFunctions = length(interpret.basisFunctions)*nodes;

% Preallocation
fisherDetMatrix = zeros(length(SNR),numRealisations,length(measurements));
mseMatrix = zeros(length(SNR),numRealisations,length(measurements),nodes);
estimate = zeros(length(SNR),numRealisations,length(measurements),numFunctions, nodes);
corrDerMatrix = zeros(length(SNR),numRealisations,1001,nodes);

% Sampling parameters from an uniform distribution interval [0,1]
parameters = ones(3,1)*[40,1,3,0.5,1];

% Create topology for the gene regulatory network
sim = geneGraph(nodes);
sim = sim.repression(1,nodes, parameters(1,1), 4);
sim = sim.degradation(1, -parameters(1,5));

for n=2:nodes
    sim = sim.repression(n,n-1, parameters(n,1),4);
    sim = sim.degradation(n, -parameters(n,5));
end

% Sampling initialConditions
initialConditions = [1; 2; 3]; % Symmetry breaking
%initialConditions = abs(10*randn(1, nodes));

% Ground truth weights
groundTruth = sim.standardGroundTruth;

[derivativeSeries, timeSeries] = ...,
    sim.runRungeKutta(initialConditions, 0, 0:0.01:10);
Phi = interpret.constructDictionary(timeSeries, false);

% As time-dictionary is same, optimal sequence has to be calculated only
% once
load sequence

% Runge Kutta simulation
for i=1:length(SNR)
    for r=1:numRealisations
        
        % IMPORTANT NOTE: Here we only corrupt the derivative series!
        [corrDer, noiseStd] =  ...,
            signalCorruption(derivativeSeries, SNR(i));
        corrDerMatrix(i,r,:,:) = corrDer;
        lambda = max(0.005, noiseStd^2);
        
        parfor j=1:length(measurements)
            % Sample data points with the maximal Fisher information
            currentIdx = allIdx{j};
            
            for l=1:nodes
                [~, estimateTemp, cost, ~, penalty, ols, convergenceGamma] = ...,
                    interpret.reconstructSetIter( ...,
                    timeSeries(currentIdx,:),...,
                    corrDer(currentIdx,l),lambda, false, 30);
                estimate(i,r,j,:,l) = estimateTemp;
                mseMatrix(i,r,j,l) = ...,
                    norm(estimateTemp-groundTruth(:,l),2)/ ...,
                    norm(groundTruth(:,l),2);
                
            end
        end
        save(['checkpoints/run_', dateChar, '_', num2str(length(measurements)), '_',  ...,
            num2str(numRealisations)]);
    end
    save(['checkpoints/run_', dateChar, '_', num2str(length(measurements)), '_',  ...,
        num2str(numRealisations)]);
end


