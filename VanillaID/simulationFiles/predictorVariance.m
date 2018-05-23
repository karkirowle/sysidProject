
function estVariance = predictorVariance(estimate, dictionary, groundTruth)
% INPUTS
% estimate - B * 1 vector containing the parameters for each dictionary
% function of the linear estimator
% dictionary - B * T vector containing the evaluated values of the
% dictionary function for the sampled x's
% groundTruth - 1 * T the ground truth derivative signal

% OUTPUT
% estVariance - 1 * 1 scalar representing the variance of the estimator

% Assumption: the estimator is unbiased
% The variance of unbiased estimator is the differece squared

estimatedValues = estimate' * dictionary;
estVariance = mean((groundTruth - estimatedValues).^2);


end
