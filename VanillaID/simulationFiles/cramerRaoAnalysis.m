% Cramer-Rao analysis
% The aim of this script is to compare Fisher information with the
% Cramer-Rao bound

% Load results
load('checkpoints/run_22_May_2018_17_56_59_50_100.mat')

% Experiment selection
SNRID = 1;
realisationID = 1;
estId = 1;

fisherVector = squeeze(fisherDetMatrix(SNRID,realisationID,:));
estimateVector = squeeze(estimate(SNRID,realisationID,:,:,estId));

T = size(fisherVector,1);
estVariance = zeros(T,1);

% Calculate prediction variance for each experiment and store in vector
for i=1:T
    sliceEstimate = estimateVector(i,:)';
    estVariance(i) = predictorVariance(sliceEstimate, Phi', derivativeSeries(:,1)');
end

figure;
%plot(log10(1./(fisherVector)));
hold on;
plot(estVariance);
run('figureFormatter');


