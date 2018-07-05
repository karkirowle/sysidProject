function [w_out, cost, w_unpruned, penalty, ols, convergenceGamma, lastGamma] = tac_reconstruction(Output, Dic, lambda,MAXITER)
%%
% This is a VANILA implementation of the follwing paper
% 
% The algorithm solves the inverse problem 
%   $y = Aw + \xi$
%
% ============= Inputs =============
% y                                    : output, in the paper $y(t) =
%                                           (x(t+\delta)-x(t))/\delta$;
% Dic                                    : dictionary matrix; 
% lambda                         : the tradeoff parameter you should use, 
%                                          basically, it is proportional to
%                                          the invese of the variance, e.g. 1;
% MAXITER                    : maximum number of the iterations,  e.g. 5
% ============ Outputs =============
% w_out - weights pruned with the heuristic
% cost - final value of the cost fuction
% penalty - finel value of the penalty part of the cost function
% (regularisation part)
% ols - final value of the Least Squares part of the cost function
% (goodness of fit)
% convergenceGamma - a boolean indicator whether the Gamma values
% have been converged or not
% lastGamma - the final values of Gamma (N*1)
%
% ============= Author =============
%  Wei Pan (w.pan11@imperial.ac.uk, panweihit@gmail.com)
% ==== Additional modifications ====
%  Bence Halpern (bmh14@ic.ac.uk, szerk.animeblog@gmail.com)
% ============= Version ============
%   1.1 (June 28, 2018)
%%

% Delta criterion to determine convergence of w
delta = 1e-6;
% Delta criterion to determine convergence of Gamma
delta2 = 1e-2;

% The dictionary is a MXN (timeXstates)
[M,N]=size(Dic); 

% Initialisation/preallocation of the variables
U=ones(N, MAXITER);
Gamma=zeros(N, MAXITER);
UU=zeros(N, MAXITER);
w_estimate=zeros(N, MAXITER);
WWW=ones(N, MAXITER);

fprintf(1, 'Sparsity optimisation ...');

for iter=1:1:MAXITER
    
    fprintf('Round %d ', iter);
    cvx_begin quiet
    cvx_solver sedumi   %sdpt3
    variable W(N)
    
    % Perform initially a LASSO, then a reweighted LASSO Optimisation 
    minimize    (lambda*norm( U(:,iter).*W, 1 )+ 0.5*sum((Dic* W-Output).^2) )
    %                 subject to
    %                           W.^2-ones(101,1)<=0;
    cvx_end
    % NOTE: you can add constraints to the problem if there is some
    % knowledge in the form of weights, uncomment above code
    
    
    % 
    w_estimate(:,iter)=W;
    WWW(:,iter)=W;
    Gamma(:,iter)=U(:,iter).^-1.*abs(W);
    
    % Obtain the MAP covariance
    Dic0=lambda*eye(M)+Dic*diag(Gamma(:,iter))*Dic';
    
    % Obtain u, which determine the reweighting coefficient based
    % on sensitivity to fit error
    UU(:, iter)=diag(Dic'*(Dic0\Dic));
    U(:,iter+1)=abs(sqrt(UU(:, iter)));
    
    w_unpruned = w_estimate(:,iter);

    % Convergence heuristic: if weights are small enough, prune them
    for i=1:N
        if   w_estimate(i,iter).^2/norm(w_estimate(:,iter))^2<delta
            w_estimate(i,iter)=0;
        end
    end
    
    % Termination of loop on convergence only if iter > 5
    if (iter > 5)
       % Claim convergence when difference of last sample and mean
       % of last five is less than some delta
       meanGamma = mean(Gamma(:,iter-4:iter),2);
       lastGamma = Gamma(:,iter);
       normGamma = norm(meanGamma - lastGamma)/norm(lastGamma);
       if (normGamma < delta2)
            convergenceGamma = true;
            break;
       end
    end
    convergenceGamma = true;

end


% Output only the weights at the last iteration
w_out = w_estimate(:, iter);

% Output the final values for the cost function
cost = lambda*norm( U(:,iter).*w_out, 1 ) + 0.5*sum((Dic* w_out-Output).^2);

penalty = lambda*norm( U(:,iter).*W, 1 );

ols = 0.5*sum((Dic* w_out-Output).^2);
