
function pass = testSNR(trueSignal,corruptedSignal, trueSNR)

% NOTE: There is a 1 unit variance in the estimated true SNR difference due
% to the random nature of noise and that we don't have infinite amount of
% timeseries

estimatedSNR = snr(trueSignal,corruptedSignal-trueSignal);
pass = norm((estimatedSNR - trueSNR),2)/norm(trueSNR,2) < 0.1 ;


end