% Bayesian Identification of Gene Regulatory Networks
% Bence Halpern 2018
% Reconstruction error analysis


% Housekeeping
clc;
clear;
close all;


% Creating two gene repressing topology
nodes = 2;
initialConditions = [20,10];
sim = geneGraph(nodes);

sim = sim.degradation(1, -0.2);
sim = sim.degradation(2, -0.2);
sim = sim.repression(1,2, 0.5,4);
sim = sim.repression(2,1,0.5,4);

% Ground truth weights -> should be extracted from geneGraph (to implement)
groundTruth = zeros(18,2);
groundTruth(1,1) = -0.2;
groundTruth(2,2) = -0.2;
groundTruth(9,2) = 0.5;
groundTruth(10,1) = 0.5;
SNR = 100;
measurements = 1:50:1000;
% Interpretation parameters
lambda = 0.01;
interpret = interpretationGraph(1); % node number does not matter here


% Adding basis functions to model
interpret = interpret.addBasisFunction(@(x) x);
for i=1:4
	interpret = interpret.addBasisFunction(@(x) 1./(1+x).^i);
end
for i=1:4
	interpret = interpret.addBasisFunction(@(x) x./(1+x).^i);
end


csvPath = [pwd ,'/results/valami.csv'];
fileID = fopen(csvPath, 'w');

title = ['Reconstruction lambda,SNR,Data length,Log MSE on Weights', newline];
fprintf(fileID, title);


% Runge Kutta simulation 
[derivativeSeries, timeSeries] = sim.runRungeKutta(initialConditions, 0, 0:0.1:200);
for k=1:length(lambda)
for i=1:length(SNR)
    [corrTime, corrDer] = signalCorruption(timeSeries,derivativeSeries,SNR(i));
    disp(['Passed:', num2str(testSNR(timeSeries(:,1), corrTime(:,1), SNR(i)))]);
    for j=1:length(measurements)
        try
            % j is amount of measurements
            % NOTE: this is not random sampling, but in this case the
            % "information gain" is becoming worse
          %  dataIndices = randi(2000,1,j);
          %  disp(dataIndices)
            disp(['Working on: lambda', num2str(lambda(k)), ...,
                ' SNR:' num2str(SNR(i)), ' DataAmount:', num2str(measurements(j))]);
            disp('Before interpret')
        [~, estimate1, ~] = interpret.reconstruct(corrTime(1:measurements(j),:), corrDer(1:measurements(j),1),lambda(k));
        [~, estimate2, ~] = interpret.reconstruct(corrTime(1:measurements(j),:), corrDer(1:measurements(j),2),lambda(k));
        estimate = [estimate1, estimate2];
        disp('Estimate done');
        % Error in weight reconstruction
        %resultsRow = {num2str(lambda(k)),  num2str(SNR(i)),  num2str(measurements(j)), num2str(mse), newline};
        %T = [T;resultsRow];
        mse = norm(estimate-groundTruth,2)/norm(groundTruth,2);
        resultsRow = [num2str(lambda(k)), ',', num2str(SNR(i)), ',', num2str(measurements(j)), ',' num2str(mse), newline];
        disp(resultsRow)
        fprintf(fileID, resultsRow);
        stdMetric(j) = std(corrTime(1:measurements(j),1));
        % TODO: Error in signal reconstruction
        catch
            disp('Failed')
        %resultsRow = {num2str(lambda(k)),  num2str(SNR(i)),  num2str(measurements(j)), 'Ill-conditioned', newline};
        %T = [T;resultsRow];
        resultsRow = [num2str(lambda(k)), ',', num2str(SNR(i)), ',', num2str(measurements(j)), ',', 'Ill-conditioned',newline];
        fprintf(fileID, resultsRow);
        end
       %disp('Round done');
    end
  
end

%fclose(fileID);
% figure(k);
% [X,Y] = meshgrid(noise,measurements);
% surf(X,Y, log(mseMatrix.'));
% xlabel('Noise parameter')
% ylabel('Number of measurements taken');
% title(['Lambda', num2str(lambda(k))])
% saveas(gcf,['Lambda', num2str(lambda(k)) '.png'])

end
