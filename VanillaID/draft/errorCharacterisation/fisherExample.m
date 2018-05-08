% Example Fisher information code

% Fisher = ecmmvnrfish(Data,Design,Covariance,Method,MatrixFormat,CovarFormat);

% Cholesky decomposition of data covariance - because it is positive
% definite -> here upper triangular and it's Hermitian L = L * L*


% Is design something you regress against??

% A = CholCovar' \ Design{k)
% which means find Design{k} * A = CholCovar'

% Mean covariance of A
% TestMatrix = TestMatrix + A'*A; 
% TestMatrix = (1.0/Count) .* TestMatrix;
% 	Fisher(1:NumParams,1:NumParams) = TestMatrix;
    
% In the notes it is consistently given as XXT/sigma^2
% The covar is an input matrix

% Housekeeping
clc;
clear;
close all;

% Covariance
data = [1:10; 1:10; 1:10];
data = data + normrnd(0,1,3,10);


% Numseries times numseries
covar = cov(data'); % numseries x numseriues for residuals?
% Numseries times numseries upper triangular
[CholCovar, p] = chol(covar);

% Numsamples by Numparams
Design = ones(10,3);
A = U' \ Design(1,:); % ez otrombaság