% Identification test function - Bence Halpern

% Housekeeping
clc;
clear;
close all;

% Tolerance parameter for unit tests
tolerance = 0.1;
testNumber = 5;

% Generate time series
[timeSeries, derivativeSeries, coefficients] = ..., 
    twoGeneNetworkGenerator(testNumber, 100, 0, 1);

% Constructing the interpreation graph
interpret = interpretationGraph(2);
interpret = interpret.addBasisFunction(@(x) x);
interpret = interpret.addBasisFunction(@(x) 1./(1+x));
interpret = interpret.addBasisFunction(@(x) 1./(1+x).^2);
interpret = interpret.addBasisFunction(@(x) 1./(1+x).^3);
interpret = interpret.addBasisFunction(@(x) 1./(1+x).^4);
interpret = interpret.addBasisFunction(@(x) x./(1+x).^1);
interpret = interpret.addBasisFunction(@(x) x./(1+x).^2);
interpret = interpret.addBasisFunction(@(x) x./(1+x).^3);
interpret = interpret.addBasisFunction(@(x) x./(1+x).^4);

for i=1:testNumber
    tempTime = reshape(timeSeries(i,1:(end-1),1:2), [(size(timeSeries,2) -1) ,2]);
    tempDeriv = derivativeSeries(i,:,1).';
    tempDeriv2 = derivativeSeries(i,:,2).';

    % Identification
    try 
    [interpret, estimate1, ~] = interpret.reconstruct(tempTime, ...,
        tempDeriv,0.1);
    [interpret, estimate2, ~] = interpret.reconstruct(tempTime, ...,
        tempDeriv2,0.1);
    estimateMatrix(i,:,:) = [estimate1, estimate2];
    
    % Tests
    test(1,i) = (estimate1(1) - coefficients(1,i)) < tolerance;
    test(2,i) = (estimate1(10) - coefficients(3,i)) < tolerance;
    test(3,i) = (estimate2(2) - coefficients(2,i)) < tolerance;
    test(4,i) = (estimate2(9) - coefficients(4,i)) < tolerance;
    catch
        test(1:4,i) = -1;
    end
end

disp("Ratio of passed tests:");
disp(sum(sum(test))/(4*testNumber));
