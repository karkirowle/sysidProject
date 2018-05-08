% Windowing type algorithm with simple repressilator
% Housekeeping
clc;
clear;
close all;

% Setting the seed
rng('default');

% Creating eight gene repressing topology
nodes = 3;
sim = geneGraph(nodes);
numRealisations = 1;

measurements = 1:1000;
lambda = 5* 10^(-3);
SNR = 0;

interpret = interpretationGraph(1); % node number does not matter here

% Adding basis functions to model
interpret = interpret.addBasisFunction(@(x) x);
for i=1:4
    interpret = interpret.addBasisFunction(@(x) 1./(1+x.^i));
end
for i=1:4
    interpret = interpret.addBasisFunction(@(x) (x.^i)./(1+x.^i));
end

numFunctions = length(interpret.basisFunctions)*nodes;

csvPath = [pwd ,'/results/repressilator3.csv'];
fileID = fopen(csvPath, 'w');

titleRow = ['Reconstruction lambda,SNR,Data length,RNMSE1,RNMSE2,Sparsity', newline];
fprintf(fileID, titleRow);

% Preallocation
mse = zeros(numRealisations,length(measurements));
mse1 = zeros(numRealisations,length(measurements));
costMatrix = zeros(numRealisations,length(measurements),nodes);
costPenalty = zeros(numRealisations,length(measurements),nodes);
costOLS = zeros(numRealisations,length(measurements),nodes);
sparsity = zeros(numRealisations,length(measurements),nodes);

estimate = zeros(numFunctions, nodes);

% Runge Kutta simulation
for r=1:numRealisations
    % Sampling parameters from an uniform distribution interval [0,1]
    %parameters = (0.8 + 0.4*rand(nodes,1)) * [40, 1, 3, 0.5, 1];
    parameters = ones(3,1)*[40,1,3,0.5,1];
    
    % Create topology for the gene regulatory network
    sim = geneGraph(nodes);
    sim = sim.repression(1,nodes, parameters(1,1), 4);
    sim = sim.degradation(1, -parameters(1,5));
    
    for i=2:nodes
        sim = sim.repression(i,i-1, parameters(i,1),4);
        sim = sim.degradation(i, -parameters(i,5));
    end
    
    % Sampling initialConditions
    initialConditions = abs(10*randn(1, nodes));
    
    % Ground truth weights
    groundTruth = sim.standardGroundTruth;
    [derivativeSeries, timeSeries] = ...,
        sim.runRungeKutta(initialConditions, 0, 0:0.01:10);
    figure;
    plot(timeSeries);
    for k=1:length(lambda)
        for i=1:length(SNR)
            corrTime = timeSeries(:,1:nodes);
            corrDer = derivativeSeries(:,1:nodes);
            for j=1:length(measurements)
                try
                    disp(['Working on: lambda', num2str(lambda(k)), ...,
                        ' SNR:' num2str(SNR(i)), ' DataAmount:', num2str(measurements(j))]);
                    for l=1:nodes
                        [~, estimateTemp, cost, ~, penalty, ols] = ...,
                            interpret.reconstructUnpruned( ...,
                            corrTime(1:measurements(j),:),...,
                            corrDer(1:measurements(j),l),lambda(k), false);
                        estimate(:,l) = estimateTemp;
                        costMatrix(r,j,l) = cost;
                        costPenalty(r,j,l) = penalty;
                        costOLS(r,j,l) = ols;
                        sparsity(r,j,l) = length(find(estimateTemp));
                    end
                    mse(r,j) = norm(estimate-groundTruth,2)/norm(groundTruth,2);
                    mse1(r,j) = norm(estimate-groundTruth,1)/norm(groundTruth,1);
                    disp(['RMNSE: ', num2str(mse(r,j))]);
                    resultsRow = [num2str(lambda(k)), ',', ...,
                        num2str(SNR(i)), ',', num2str(measurements(j)),...,
                        ',' num2str(mse(r,j)), ',' num2str(mse1(r,j)), ',' ...,
                        num2str(sparsity(r,j)), newline];
                    disp(resultsRow);
                    fprintf(fileID, resultsRow);
                catch
                    disp('Failed')
                    resultsRow = [num2str(lambda(k)), ',', ...,
                        num2str(SNR(i)), ',', num2str(measurements(j)), ...,
                        ',', 'Ill-conditioned',newline];
                    fprintf(fileID, resultsRow);
                end
            end
        end
    end
end


