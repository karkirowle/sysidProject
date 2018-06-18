% Gives back the RNMSE of the oracle estimator
% Phi - dictionary matrix k * l
% SNR - signal to noise ratio
% groundTruth - l*1 groundtruth weight vector
% derivativeSeries - k * 1 vector
% measurements - number of measurements (50)

function oracleRNMSE = oracleCalculator(Phi, SNR, groundTruth, ...,
    derivativeSeries, measurements, allCorrDer)

% Realisation parameter
realisation = size(allCorrDer,1);
% Obtain subset from ground truth
trueSubset = find(groundTruth ~= 0);

% Calculate idx
[~, noiseStd] = signalCorruption(derivativeSeries, SNR);
[~, idx] = maxRegDetDictionaryBatch(Phi', noiseStd.^2, ...,
    measurements);

% Construct subset dictionary
oracleRNMSEReal = zeros(length(idx), realisation);
for i=1:length(idx)
    subsetPhi = Phi(idx(1:i),trueSubset);
    
    % Calculate least squares
    % l * 1
    for r=1:realisation
    oracleEstimate = pinv(subsetPhi' * subsetPhi) * subsetPhi' * allCorrDer(r,idx(1:i))';
    
    % Calculate RNMSE
    oracleEstimateVector = zeros(27,1);
    oracleEstimateVector(trueSubset,:) = oracleEstimate;
    oracleRNMSEReal(i,r) = norm(oracleEstimateVector-groundTruth,2)./norm(groundTruth,2);
    end
end

oracleRNMSE = mean(oracleRNMSEReal,2);