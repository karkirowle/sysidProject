function [corruptedTimeSeries, corruptedDerivative] =  ...,
    signalCorruption(timeSeries,derivativeSeries, SNR)
    corruptedTimeSeries = awgn(timeSeries,SNR,'measured');
    corruptedDerivative = awgn(derivativeSeries,SNR,'measured');
end