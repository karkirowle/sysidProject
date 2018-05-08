% RMSE based measurement adding algorithm to investigate data variance
% For 8th differential equaiton

% Housekeeping
clc;
clear;
close all;


% Setting the seed
rng('default');

% Creating eight gene repressing topology
nodes = 8;
sim = geneGraph(nodes);


% Sampling parameters from an uniform distribution interval [0,1]

parameters = (0.8 + 0.4*rand(8,1)) * [40, 1, 3, 0.5, 1];

% Create topology for the gene regulatory network
sim = sim.repression(1,8, parameters(1,1), 3);
sim = sim.degradation(1, -parameters(1,5));

for i=2:8
    sim = sim.repression(i,i-1, parameters(i,1),3);
    sim = sim.degradation(i, -parameters(i,5));
end

% Sampling initialConditions
initialConditions = rand(1,8);

% Ground truth weights -> should be extracted from geneGraph (to implement)
groundTruth = zeros(72,8);
for i=1:8
    groundTruth(i,i) = -parameters(i,5);
end
groundTruth(32,1) = parameters(1,1);
for i=2:8
    groundTruth(23+i,i) = parameters(i,1);
end


measurements = 1:100;

lambda = 0.1;
interpret = interpretationGraph(1); % node number does not matter here


% Adding basis functions to model
interpret = interpret.addBasisFunction(@(x) x);
for i=1:4
    interpret = interpret.addBasisFunction(@(x) 1./(1+x).^i);
end
for i=1:4
    interpret = interpret.addBasisFunction(@(x) x./(1+x).^i);
end

% Runge Kutta simulation
[derivativeSeries, timeSeries] = sim.runRungeKutta(initialConditions, 0, 0:0.001:200);

% Let's take 5 subsets of the data
% Quantify the mean data variance (timeSeries std)

numSets = 10;
measurementSet = cell(numSets,1);
for i=1:numSets
    measurementSet{i}.id = 1:(100*i);
    measurementSet{i}.var = mean(var(timeSeries(measurementSet{i}.id,:)));
end

costVector = zeros(numSets,length(measurementSet{numSets}.id));
penaltyVector = zeros(numSets,length(measurementSet{numSets}.id));
olsVector = zeros(numSets,length(measurementSet{numSets}.id));
mseError = zeros(numSets,length(measurementSet{numSets}.id));
sparsityCount = zeros(numSets,length(measurementSet{numSets}.id));


for j=1:numSets
    % Randomly sample 2
    allMeasurements = measurementSet{j}.id;
    sampledIndices = [];
    id = randi(length(allMeasurements),1,2);
    sampledIndices = [sampledIndices id];
    allMeasurements(id) = [];
    
    % First reconstruction of 8
    [~, estimate, cost, estimateUnpruned, penalty, ols] = ...,
        interpret.reconstructUnpruned(timeSeries(sampledIndices,:), ...,
        derivativeSeries(sampledIndices,8),lambda);
    costVector(j,1) = cost;
    costIndex = 1;
    while (~isempty(allMeasurements))
        disp(length(allMeasurements));
        clear MSE;
        MSE = Inf*ones(1,length(allMeasurements));
        for i=1:length(allMeasurements)
            Phi = interpret.constructDictionary( ...,
                timeSeries(allMeasurements(i),:));
            MSE(i) = immse(Phi*estimate,derivativeSeries(allMeasurements(i),8));
        end
        [~, id] = max(MSE);
        sampledIndices = [sampledIndices id];
        allMeasurements(id) = [];
        [~, estimate, cost, estimateUnpruned, penalty, ols] = ...,
            interpret.reconstructUnpruned(timeSeries(sampledIndices,:), ...,
            derivativeSeries(sampledIndices,8),lambda);
        costIndex = costIndex + 1;
        costVector(j,costIndex) = cost;
        penaltyVector(j,costIndex) = penalty;
        olsVector(j,costIndex) = ols;
        mseError(j,costIndex) = norm(estimate-groundTruth,2)/norm(groundTruth,2);
        sparsityCount(j,costIndex) = length(find(estimate));
        % Time series is the basis substitued  2x72 72x1
        % estimateTemp*timeSeries(sampledIndices,:)
    end
end

figure;
plot(costVector');
title('Optimisation end costs of different measurement subsets')
xlabel('number of measurements (with maximal RMSE')
ylabel('cost');
legend(['\sigma_d = ', num2str(measurementSet{1}.var)], ...,
    ['\sigma_d = ', num2str(measurementSet{2}.var)], ...,
     ['\sigma_d = ', num2str(measurementSet{3}.var)], ...,
      ['\sigma_d = ', num2str(measurementSet{4}.var)], ...,
       ['\sigma_d = ', num2str(measurementSet{5}.var)]);
   
    
