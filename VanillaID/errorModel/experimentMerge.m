% Bayesian Identification of Gene Regulatory Networks
% Bence Halpern 2018

% The following script loads SNR results from a list of mat files and
% combines them

% Housekeeping
clc;
clear;
close all;

% List of mat files to combine together
fileNames = {'checkpoints/run_13_Jun_2018_10_41_42_50_200.mat', ...,
    'checkpoints/run_12_Jun_2018_19_54_28_50_200.mat'};
numFiles = size(fileNames,2);
savedFileName = 'MAPdata';


% Parameters for sanity check
% Realisations
r = 200;
% Measurements
m = 50;
% States
s = 3;

allSNR = [];
estimateTemp = [];
allFisher = [];
allEstimate = [];
allCorrDer = [];

for i=1:numFiles
    % Load mat file
    load(fileNames{i},'-mat','SNR','Phi', ...,
        'mseMatrix','fisherDetMatrix','estimate','groundTruth','derivativeSeries', ...,
        'corrDerMatrix');
    
    % EXCEPTION RULES FOR THESE EXPERIMENTS!!!
    %     if (i==1)
    %         SNR = SNR(1:end-2);
    %         mseMatrix = mseMatrix(1:end-2,:,:,:);
    %         fisherDetMatrix = fisherDetMatrix(1:end-2,:,:,:);
    %         estimate = estimate(1:end-2,:,:,:,:);
    %         corrDerMatrix = corrDerMatrix(1:end-2,:,:,:);
    %
    %     end
    if (i==2)
        SNR = SNR(2:end-2);
        mseMatrix = mseMatrix(2:end-2,:,:,:);
        fisherDetMatrix = fisherDetMatrix(2:end-2,:,:,:);
        estimate = estimate(2:end-2,:,:,:,:);
        corrDerMatrix = corrDerMatrix(2:end-2,:,:,:);
    end
    
    % Extract SNR
    allSNR = cat(2,SNR, allSNR);
    assert(length(SNR) == size(mseMatrix,1));
    % Order is only preserved if its the same order as SNR
    estimateTemp = cat(1,mseMatrix,estimateTemp);
    allFisher = cat(1,fisherDetMatrix,allFisher);
    allEstimate = cat(1,estimate,allEstimate);
    allCorrDer = cat(1,corrDerMatrix,allCorrDer);
    
    
end

% Do a descending ordering SNR
[totalSNR, sortIdx] = sort(allSNR, 'ascend');

% Construct totalMse matrix but first check everything is as expected
totalMse = estimateTemp(sortIdx,:,:,:);
allFisher = allFisher(sortIdx,:,:,:);
allEstimate = allEstimate(sortIdx,:,:,:,:);
allCorrDer = allCorrDer(sortIdx,:,:,:);

% Plot result to check everything is right
figure;
plot(1:m, squeeze(mean(totalMse(:,:,:,1),2)), 'LineWidth', 1.5);
legend({num2str(totalSNR')});

% Save results
clearvars -except totalSNR totalMse allFisher allCorrDer groundTruth derivativeSeries allEstimate savedFileName Phi
save(savedFileName);
