
% Bayesian Identification of Gene Regulatory Networks

% Bence Halpern 2018

% ------------------------- DESCRIPTION -----------------------------------
% Type II maximum a posteriori covariance based optimal data sampling based
% repressilator model for testing identification lower bounds

% Housekeeping
clc;
clear;
close all;

% ------------------- PARAMETERS: CHANGE WISELY! --------------------------
measurements = 1:50;
SNR = [0.1, 1, 10, 100];
numExperiments = length(SNR);
numRealisations = 200;

% Server parameters
% clusterNumber = 2;
% parpool('local',clusterNumber)

% -------------------------------------------------------------------------

% Fetching date for filename
filenameDate = datetime('now','Format','default');

% Replacing some characters to underscores so that every file name is legible
dateChar = char(filenameDate);
dateChar(dateChar == ' ') = '_';
dateChar(dateChar == '-') = '_';
dateChar(dateChar == ':') = '_';


% Setting the seed to ensure reproducibility of experimental results
rng('default');

% Creating three gene repressilator topology
nodes = 3;
sim = geneGraph(nodes);

% Don't change the one argument here!!
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
covDetMatrix = zeros(numExperiments,numRealisations,length(measurements));
mseMatrix = zeros(numExperiments,numRealisations,length(measurements),nodes);
estimate = zeros(numExperiments,numRealisations,length(measurements),numFunctions, nodes);
corrDerMatrix = zeros(numExperiments,numRealisations,1001,nodes);

% Runge Kutta simulation
for i=1:numExperiments
    estimateRTemp = zeros(numRealisations,length(measurements),nodes, numFunctions);
    mseRTemp = zeros(numRealisations,length(measurements),nodes);
    
    for r=1:numRealisations
        % Sampling parameters from an uniform cdistribution interval [0,1]
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
        
        % Ground truth weights
        groundTruth = sim.standardGroundTruth;
        
        [derivativeSeries, timeSeries] = ...,
            sim.runRungeKutta(initialConditions, 0, 0:0.01:10);
        
        % IMPORTANT NOTE: Here we only corrupt the derivative series!
        [corrDer, noiseStd] =  ...,
            signalCorruption(derivativeSeries, SNR(i));
        corrDerMatrix(i,r,:,:) = corrDer;
        lambda = max(0.005, noiseStd^2);
        
        Phi = interpret.constructDictionary(timeSeries, false);
        estimateMTemp = zeros(length(measurements),nodes,numFunctions);
        mseMTemp = zeros(length(measurements),nodes);
        for j=1:length(measurements)
            
            % Edge case at start
            if (j == 1)
                PhiP = [];
                currentIdx = [];
                Gamma = ones(1,numFunctions);
            end
            
            % Sample data points with the maximal Fisher information
            [gammaInfo, PhiP, currentIdx] = maxGammaDictionary(Phi', lambda, ...,
                PhiP, currentIdx, diag(Gamma));
            
            estimateBTemp = zeros(nodes,numFunctions);
            mseBTemp = zeros(1,nodes);
            for l=1:nodes
                [~, estimateTemp, cost, ~, penalty, ols, convergenceGamma,Gamma] = ...,
                    interpret.reconstructGamma( ...,
                    timeSeries(currentIdx,:),...,
                    corrDer(currentIdx,l),lambda, false, 30);
                %estimate(i,r,j,:,l) = estimateTemp;
                estimateBTemp(l,:) = estimateTemp;
                
%                 mseMatrix(i,r,j,l) = ...,
%                     norm(estimateTemp-groundTruth(:,l),2)/ ...,
%                     norm(groundTruth(:,l),2);
                mseBTemp = norm(estimateTemp-groundTruth(:,l),2)/ ...,
                    norm(groundTruth(:,l),2);
            end
            estimateMTemp(j,:,:) = estimateBTemp;
            mseMTemp(j,:) = mseBTemp;
            
        end
        estimateRTemp(r,:,:,:) = estimateMTemp;
        mseRTemp(r,:,:) = mseMTemp;
    end
    estimate(i,:,:,:,:) = estimateRTemp;
    mseMatrix(i,:,:,:) = mseRTemp;
    save(['checkpoints/run_', dateChar, '_', num2str(length(measurements)), '_',  ...,
        num2str(numRealisations)]);
end


