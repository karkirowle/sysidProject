% Represisng motif test function
% Housekeeping
clc;
clear;
close all;

% Testnumber parameter
testNumber = 10;
lambda = 0.1;
% Generate time series
[timeSeries, derivativeSeries, coefficients] = ..., 
    twoGeneNetworkGenerator(testNumber, 100, 0, 1);

% Constructing the interpreation graph for the repressing motif
interpretRepress = interpretationGraph(1);
interpretRepress = interpretRepress.addBasisFunction(@(x) x);
interpretRepress = interpretRepress.addBasisFunction(@(x) 1./(1+x));
interpretRepress = interpretRepress.addBasisFunction(@(x) 1./(1+x).^2);
interpretRepress = interpretRepress.addBasisFunction(@(x) 1./(1+x).^3);
interpretRepress = interpretRepress.addBasisFunction(@(x) 1./(1+x).^4);

% Constructin the interpretation graph for the activating motif
interpretActivate = interpretationGraph(1);
interpretActivate = interpretActivate.addBasisFunction(@(x) x);
interpretActivate = interpretActivate.addBasisFunction(@(x) x./(1+x).^1);
interpretActivate = interpretActivate.addBasisFunction(@(x) x./(1+x).^2);
interpretActivate = interpretActivate.addBasisFunction(@(x) x./(1+x).^3);
interpretActivate = interpretActivate.addBasisFunction(@(x) x./(1+x).^4);

testResults = zeros(testNumber,1);

for i=1:testNumber
    % NOTE: This is taking only the first measurement set! 
    tempTime = reshape(timeSeries(i,1:(end-1),1), [(size(timeSeries,2) -1),1]);
    tempDeriv = derivativeSeries(i,:,1).';

    % Identification

    [interpretRepress, ~, lossRep(i)] = interpretRepress.reconstruct(tempTime, ...,
        tempDeriv,lambda);
 

    [interpretActivate, ~, lossAct(i)] = interpretActivate.reconstruct(tempTime, ...,
        tempDeriv,lambda);
 
    testResults(i) = lossAct(i) < lossRep(i);
end

disp("Ratio of passed tests:");
disp(sum(testResults)/length(testResults));