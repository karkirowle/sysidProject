% Adds measurement noise to the time series and the derivative series
function [corruptedTimeSeries, corruptedDerivative, noiseStd] =  ...,
    signalCorruption(timeSeries,derivativeSeries, SNR)

% White gaussian noise added to time series
corruptedTimeSeries = awgn(timeSeries,SNR,'measured');

% White gaussian noise added to derivative series
corruptedDerivative = awgn(derivativeSeries,SNR,'measured');

% Signal to noise ratio is measured std / noise std
sigPower = sum(abs(timeSeries(:)).^2)/length(timeSeries(:));
noisePower = sigPower/SNR;
noiseStd = sqrt(noisePower);
end