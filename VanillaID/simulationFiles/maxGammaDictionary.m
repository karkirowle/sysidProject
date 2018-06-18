function [gammaInfos, idx] = maxGammaDictionary(data, lambda, ...,
    numGamma, previousPhi)
% INPUT
% data - P * T matrix which contains the values of the dictionary at T time
% lambda - constant for calculating Fisher information points
% previousPhi - dictionary from the previous iteration of the algorithm, if
% empty will assume this is the first iteration

% OUTPUTS
% fisherInfo - T * 1 matrix which constains Fisher Information value witch each
% added points
% idx - T*1 list of indices which sequentially maximise Fisher information

% Check variance is positive
assert(lambda >= 0);

% Useful constants
P = size(data,1);
T = size(data,2);

% Preallocation
idx = zeros(numGamma,1);
fisherInfos = zeros(numGamma,1);

% Check if previousPhi is empty
if (isempty(previousPhi))
    
else
    
end



for i=1:numGamma
    
    % Treatment of edge case where there aren't any selected idx yet
    if (i==1)
        selectedIdx = [];
    else
        selectedIdx = idx(1:i-1);
    end
    
    % Select set of indices above which we desire to perform max search
    indices = (1:T)';
    indices(selectedIdx) = [];
    
    
    % Calculate fisher information for each combination
    gammaSearch = zeros(length(indices),1);
    for j=1:length(indices)
        X = data(:,[selectedIdx;indices(j)]);
        % The determinant is taken here to obtain a scalar (D-optimal)
        gammaSearch(j) = det(Gamma\(lambda*eye(P)) + X * X');
    end
    
    % Find maximal Fisher information in the above combination vector
    [maxFisher,maxId] = max(gammaSearch);
    fisherInfos(i) = maxFisher;
    idx(i) = indices(maxId);
end

end