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
SNR = [50,40,30,20,10];
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

% Runge Kutta simulation
for i=1:length(SNR)
    for r=1:numRealisations
        % Sampling parameters from an uniform distribution interval [0,1]
        %parameters = (0.8 + 0.4*rand(nodes,1)) * [40, 1, 3, 0.5, 1];
        parameters = ones(3,1)*[40,1,3,0.5,1];
        
        % For the first differential equation now one is sampled randomly
        %parameters(1,5) = (0.8 + 0.4*rand) * 40;
        
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
        % IMPORTANT NOTE: Here we only corrupt the derivative series!
        [corrDer, noiseStd] =  ...,
            signalCorruption(derivativeSeries, SNR(i));
        corrDerMatrix(i,r,:,:) = corrDer;
        lambda = max(0.005, noiseStd^2);
        
        % Calculate batch ordering of maximal Fisher information
        Phi = interpret.constructDictionary(timeSeries, false);
        [fisherInfos, idx] = maxSupervisedFisherDictionary(Phi', 10^(-5), ...,
            corrDer', length(measurements));
        fisherDetMatrix(i,r,:) = fisherInfos;
        parfor j=1:length(measurements)
            % Sample data points with the maximal Fisher information
            currentIdx = idx(1:j);
            
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


