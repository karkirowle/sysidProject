% Reproduction based on he error curves of the paper:
% Identifying Biochemical Reaction Networks From Heterogeneous Datasets
% Data analysis done by Bence Halpern

% He takes more realisations and averages them - 50 namely - CHANGED
% Initial conditions are also multiplied by times 10 - CHANGED 
% He used Euler method - we will use RungeKutta, because Euler is SLOW

% Housekeeping
clc;
clear;
close all;

% Setting the seed
rng('default');

% Creating eight gene repressing topology
nodes = 8;


SNR = 1000;
measurements = 1:100;

% Interpretation parameters
lambda = 5 * 10^(-3);
numRealisations = 50;
interpret = interpretationGraph(1); % node number does not matter here


% Adding basis functions to model
interpret = interpret.addBasisFunction(@(x) x);
interpret = interpret.addBasisFunction(@(x) 1./(1+x.^3));
interpret = interpret.addBasisFunction(@(x) (x.^3)./(1+x.^3));



csvPath = [pwd ,'/results/heteroReproductionWeiCodeBased.csv'];
fileID = fopen(csvPath, 'w');

titleRow = ['Reconstruction lambda,SNR,Data length,MNSE', newline];
fprintf(fileID, titleRow);

% Preallocation
mse = zeros(numRealisations,length(measurements));
estimate = zeros(25, nodes);

% Runge Kutta simulation 
for r=1:numRealisations
    % Sampling parameters from an uniform distribution interval [0,1]
    
    parameters = (0.8 + 0.4*rand(8,1)) * [40, 1, 3, 0.5, 1];
    
    % Create topology for the gene regulatory network
    sim = geneGraph(nodes+1); %
    sim = sim.repression(1,8, parameters(1,1), 3);
    sim = sim.degradation(1, -parameters(1,5));
    sim = sim.bias(1, 9, parameters(1,4));
    
    for i=2:8
        sim = sim.repression(i,i-1, parameters(i,1),3);
        sim = sim.degradation(i, -parameters(i,5));
        sim = sim.bias(i, 9, parameters(i,4));

    end
    
    % Sampling initialConditions
    initialConditions = [abs(10*randn(1, nodes)), 1]; % bias in cond
    
    % Ground truth weights -> should be extracted from geneGraph (to implement)
    groundTruth = zeros(25,8);
    for i=1:8
        groundTruth(i,i) = -parameters(i,5);
    end
    groundTruth(16,1) = parameters(1,1);
    for i=2:8
        groundTruth(7+i,i) = parameters(i,1);
    end
    groundTruth(25,:) = parameters(:,4)';
    [derivativeSeries, timeSeries] = sim.runRungeKutta(initialConditions, 0, 0:0.1:10);
    for k=1:length(lambda)
        for i=1:length(SNR)
            %[corrTime, corrDer] = signalCorruption(timeSeries,derivativeSeries,SNR(i));
            corrTime = timeSeries(:,1:8);
            corrDer = derivativeSeries(:,1:8);
            for j=[1, 10:10:length(measurements)]
                try
                    disp(['Working on: lambda', num2str(lambda(k)), ...,
                        ' SNR:' num2str(SNR(i)), ' DataAmount:', num2str(measurements(j))]);
                    for l=1:nodes
                        [~, estimateTemp, ~, ~, ~, ~] = ...,
                            interpret.reconstructUnpruned( ...,
                            corrTime(1:measurements(j),:),...,
                            corrDer(1:measurements(j),l),lambda(k), true);
                        estimate(:,l) = estimateTemp;
                    end
                    mse(r,j) = norm(estimate-groundTruth,2)/norm(groundTruth,2);
                    disp(['MNSE: ', num2str(mse(r,j))]);
                    resultsRow = [num2str(lambda(k)), ',', ...,
                        num2str(SNR(i)), ',', num2str(measurements(j)),...,
                        ',' num2str(mse(r,j)), newline];
                    disp(resultsRow)
                    fprintf(fileID, resultsRow);
                catch
                    disp('Failed')
                    resultsRow = [num2str(lambda(k)), ',', ...,
                        num2str(SNR(i)), ',', num2str(measurements(j)), ...,
                        ',', 'Ill-conditioned',newline];
                    fprintf(fileID, resultsRow);
                end
            end
            figure;
            plot(1:100, mean(mse(r,j),1));
            title('RNMSE plot')
            xlabel('measurement amount');
            ylabel('RNMSE averaged over 50 realisations');
        end
    end
end


