% Adds measurement noise to the time series and the derivative series
function [corruptedDerivative, noiseStd] =  ...,
    signalCorruption(derivativeSeries, SNR)

% Parameters
M = size(derivativeSeries,1);
N = size(derivativeSeries,2);

% Preallocation
sigPower = zeros(1,N);

% Noise power is estimated for each channel of derivatives
for i=1:N
    sigPower(i) = sum(abs(derivativeSeries(:,i)).^2)/length(derivativeSeries(:,i));
end

% Warning if noiseStd are very different for each channel
if (max(sigPower) - min(sigPower) > 100)
   warning('Channel signal powers vary largely. SNR estimate might be inaccurate!');
end

% Signal to noise ratio is measured std / noise std
noisePower = mean(sigPower,2)./SNR;
noiseStd = sqrt(noisePower);

% Noise sampled for a WGN
corruptedDerivative = derivativeSeries + noiseStd.*randn(M,N);
end