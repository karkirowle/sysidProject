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
measurements = 1:10;
SNR = 50;
numRealisations = 1;


% -------------------------------------------------------------------------

% Fetching date for filename
filenameDate = datetime('now','Format','default');
dateChar = char(filenameDate);
dateChar(dateChar == ' ') = '_';
dateChar(dateChar == '-') = '_';
dateChar(dateChar == ':') = '_';


% Setting the seed to ensure reproducibility of experimental results
rng('default');

for b=5:10
    
    % Creating three gene repressilator topology
    nodes = 3;
    sim = geneGraph(nodes);
    
    % Don't change the one argument here
    interpret = interpretationGraph(1);
    
    % Adding basis functions to the repressilator model
    
    % Degradation function (linear)
    interpret = interpret.addBasisFunction(@(x) x);
    
    % Hill functions added up to order 4
    for i=1:b
        interpret = interpret.addBasisFunction(@(x) 1./(1+x.^i));
    end
    
    for i=1:b
        interpret = interpret.addBasisFunction(@(x) (x.^i)./(1+x.^i));
    end
    
    
    numFunctions = length(interpret.basisFunctions)*nodes;
    
    % Preallocation
    fisherDetMatrix = zeros(numRealisations,length(measurements));
    mseMatrix = zeros(numRealisations,length(measurements),nodes);
    estimate = zeros(numRealisations,length(measurements),numFunctions, nodes);
    
    % Runge Kutta simulation
    i = 1;
    for r=1:numRealisations
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
        
        % Ground truth weights dependent on number of functions
        groundTruth = zeros((2*b+1)*3,3);
        groundTruth(1,1) = -1;
        groundTruth(2,2) = -1;
        groundTruth(3,3) = -1;
        groundTruth(3*b + 3, 1) = 40;
        groundTruth(3*b + 1, 2) = 40;
        groundTruth(3*b + 2, 3) = 40;
        
        % IMPORTANT NOTE: Here we only corrupt the time series!
        [derivativeSeries, timeSeries] = ...,
            sim.runRungeKutta(initialConditions, 0, 0:0.01:10);
        [corrDer, noiseStd] =  ...,
            signalCorruption(derivativeSeries, SNR);
        lambda = max(0.005, noiseStd^2);
        
        % Calculate batch ordering of maximal Fisher information
        Phi = interpret.constructDictionary(timeSeries, false);
        [fisherInfos, idx] = maxFisherDictionaryBatch(Phi', lambda, ...,
            length(measurements));
        fisherDetMatrix(r,:) = fisherInfos';
        figure;
        plot(derivativeSeries);
        parfor j=1:length(measurements)
            % Sample data points with the maximal Fisher information
            currentIdx = idx(1:j);
            
            for l=1:nodes
                [~, estimateTemp, cost, ~, penalty, ols, convergenceGamma] = ...,
                    interpret.reconstructSetIter( ...,
                    timeSeries(currentIdx,:),...,
                    corrDer(currentIdx,l),lambda, false, 30);
                estimate(r,j,:,l) = estimateTemp;
                mseMatrix(r,j,l) = ...,
                    norm(estimateTemp-groundTruth(:,l),2)/ ...,
                    norm(groundTruth(:,l),2);
                
            end
            disp('Pakk');
            
        end
        fisherDetMatrixStruct{b} = fisherDetMatrix;
        mseMatrixStruct{b} = mseMatrix;
        estimateStruct{b} = estimate;
        save(['checkpoints/run_', dateChar, '_', num2str(length(measurements)), '_',  ...,
            num2str(numRealisations)]);
    end
    save(['checkpoints/run_', dateChar, '_', num2str(length(measurements)), '_',  ...,
        num2str(numRealisations)]);
end

