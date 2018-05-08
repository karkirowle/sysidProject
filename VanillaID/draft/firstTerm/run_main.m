
clear all;
close all;

load('rep_exp1_bence.mat')
which = 1;
noise = 0.01; % standard deviation of the noise
param = getGRNParameters;
[y, A, w_true] = GRN_dis(which,noise,param);
% [y, A, w_true] = GRN_con(which);

%w_lasso = lasso(A,y);
%w_lasso = w_lasso(:,end);
%w_ls = lscov(A,y); 
lambda = 0.01; MAXITER = 5;

% The construction of the dictionary indicates how many states it uses I
% think
w_ours =  tac_reconstruction(x_2(:,4), A, lambda,MAXITER);

%compare = [w_true, w_ours(:,end), w_lasso, w_ls]
