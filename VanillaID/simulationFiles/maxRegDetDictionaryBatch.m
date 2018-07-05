function [fisherInfos, idx] = maxRegDetDictionaryBatch(data, lambda, numFisher)
% INPUT
% data - P * T matrix which contains the values of the dictionary at T time
% lambda - constant for calculating Fisher information points
% numFisher - maximal number of points we wish to consider
% OUTPUTS
% fisherInfo - T * 1 matrix which constains Fisher Information value witch each
% added points
% idx - T*1 list of indices which sequentially maximise Fisher information


% Useful constants
P = size(data,1);
T = size(data,2);

% Preallocation
idx = zeros(numFisher,1);
fisherInfos = zeros(numFisher,1);

for i=1:numFisher
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
    fisherSearch = zeros(length(indices),1);
    for j=1:length(indices)
        X = data(:,[selectedIdx;indices(j)]);
        % The determinant is taken here to obtain a scalar (D-optimal)
        fisherSearch(j) = det(lambda*eye(P) + X * X');
    end
    
    % Find maximal Fisher information in the above combination vector
    [maxFisher,maxId] = max(fisherSearch);
    fisherInfos(i) = maxFisher;
    idx(i) = indices(maxId);
end

end