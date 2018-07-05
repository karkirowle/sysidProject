
function [F,idx] = fisherOptimalAdding(data, numFisher, lambda)
% numFisher - number of measurements to use in the Fisher information
% matrix
% data - Time*Dictionary functions matrix
% idx - samples to add
% Function based on Boyd's CVX Example

% Shorthand for time
p = size(data,1);

cvx_begin quiet
  variable lambda(p)
  maximize ( det_rootn( data'*diag(lambda)*data ))
  subject to
    sum(lambda) == numFisher;
    lambda >= 0;
    lambda <= 1;
cvx_end

lambdaD = lambda; 
[~, allIdx] = sort(lambdaD, 'descend');
idx = allIdx(1:numFisher);
F = det( data(idx,:)' * data(idx,:));

end