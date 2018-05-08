function [F, I] = maxFisherDictionary(data, lambda, previousIndices)
% INPUT
% x - P * T matrix which contains the values of the dictionary at T time
% lambda - constant for calculating Fisher information
% points
% previousIndices - previousMaxIndices (randomly sampled first index for
% first iteration
% OUTPUTS
% F - Fisher Information value
% I - data points indices with max fisher information

% TODO: optionally other optimality techniques

indices = 1:size(data,2);
indices(previousIndices) = [];

fisher = zeros(length(indices),1);
for i=1:length(indices)
    X = data(:,[previousIndices,indices(i)]);
    fisher(i) = det(X * X' ./ lambda^2);
end

[F,I1] = max(fisher);
I2 = indices(I1);
I = [previousIndices, I2];

end