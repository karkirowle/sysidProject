% The following script loads SNR results from a list of mat files and
% combines them

% Housekeeping
clc;
clear;
close all;

% List of mat files to combine together
fileNames = {'checkpoints/run_14_Jun_2018_07_11_15_50_200.mat'};
numFiles = size(fileNames,2);
savedFileName = 'Randomdata';


% Parameters for sanity check
% Realisations
r = 200;
% Measurements
m = 50;
% States
s = 3;

noiseSNR = [];
estimateTemp = [];

noiseEstimate = [];
noiseCorrDer = [];

for i=1:numFiles
    % Load mat file
    load(fileNames{i},'-mat','SNR','Phi', ...,
        'mseMatrix','estimate','groundTruth','derivativeSeries', ...,
        'corrDerMatrix');
    

    
    % Extract SNR
    noiseSNR = cat(2,SNR, noiseSNR);
    assert(length(SNR) == size(mseMatrix,1));
    % Order is only preserved if its the same order as SNR
    estimateTemp = cat(1,mseMatrix,estimateTemp);
    noiseEstimate = cat(1,estimate,noiseEstimate);
    noiseCorrDer = cat(1,corrDerMatrix,noiseCorrDer);
    
    
end

% Do a descending ordering SNR
[totalNoiseSNR, sortIdx] = sort(noiseSNR, 'ascend');

% Construct totalMse matrix but first check everything is as expected
totalNoiseMse = estimateTemp(sortIdx,:,:,:);
noiseEstimate = noiseEstimate(sortIdx,:,:,:,:);
noiseCorrDer = noiseCorrDer(sortIdx,:,:,:);

% Plot result to check everything is right
figure;
plot(1:m, squeeze(mean(totalNoiseMse(:,:,:,1),2)), 'LineWidth', 1.5);
legend({num2str(totalNoiseSNR')});

% Save results
clearvars -except totalNoiseSNR totalNoiseMse noiseCorrDer groundTruth derivativeSeries noiseEstimate savedFileName Phi
save(savedFileName);
