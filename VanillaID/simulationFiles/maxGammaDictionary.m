function [gammaInfo, newPhi, idx] = maxGammaDictionary(data, lambda, ...,
    previousPhi, previousIdx, Gamma)
% INPUT
% data - P * T matrix which contains the values of the dictionary at T time
% lambda - constant for calculating Fisher information points
% previousPhi - L*T where L < P dictionary from the previous iteration of the algorithm, if
% empty will assume this is the first iteration
% previousIdx - L*1
% Gamma - P * P
% OUTPUTS
% gammaInfo - T * 1 matrix which constains Fisher Information value witch each
% added points
% newPhi - (L+1)*T new dictionary
% idx - T*1 list of indices which sequentially maximise Fisher information

% Check variance is positive
assert(lambda >= 0);
% Check that the selectedIdx is the same damension
assert(size(previousPhi,2) == size(previousIdx,1),['Previous dcitionary size has ', ...,
    num2str(size(previousPhi,2)), ' rows while the respective ID vector have only ', ...,
    num2str(size(previousIdx,1))]);

% Useful constants
P = size(data,1);
T = size(data,2);
L = size(previousIdx,1);



% Check if previousPhi is empty
if (L == 0)
    selectedIdx = [];
else
    selectedIdx = previousIdx;
end

% Select set of indices included in the maximum search
indices = (1:T)';
indices(selectedIdx) = [];

% Calculate the informativeness measure, first preallocate
gammaSearch = zeros(length(indices),1);

for j=1:length(indices)
    
    % Concatenate with previous result
    X = data(:,[selectedIdx;indices(j)]);
    
    % The determinant is taken here to map the information to a scalar (D-optimal)
    gammaSearch(j) = det(pinv(Gamma)*(lambda*eye(P)) + X * X');
    
end

% Find maximal Fisher information in the above combination vector
[maxGamma,maxId] = max(gammaSearch);
gammaInfo= maxGamma;

% Add the newfound id
idx = [previousIdx; indices(maxId)];

% Return the new dictionary
newPhi = data(:,[selectedIdx;indices(maxId)]);
end