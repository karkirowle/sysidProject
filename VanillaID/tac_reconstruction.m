function [w_out, cost, w_unpruned, penalty, ols, convergenceGamma] = tac_reconstruction(Output, Dic, lambda,MAXITER)
%%
% This is a VANILA implementation of the follwing paper
% 
% The algorithm solves the inverse problem 
%   $y = Aw + \xi$
%
% ============= Inputs ===============
% y                                    : output, in the paper $y(t) =
%                                           (x(t+\delta)-x(t))/\delta$;
% Dic                                    : dictionary matrix; 
% lambda                         : the tradeoff parameter you should use, 
%                                          basically, it is proportional to
%                                          the invese of the variance, e.g. 1;
% MAXITER                    : maximum number of the iterations, e.g. 5

% ============= Reference =============
% W. Pan, Y. Yuan, J. Goncalves, and G.-B. Stan, 
% A Sparse Bayesian Approach to the Iden- tification of Nonlinear State-Space Systems,
% IEEE Transaction on Automatic Control, 2015 (to appear). arXiv:1408.3549
% http://arxiv.org/abs/1408.3549
%
% ============= Author =============
%  Wei Pan (w.pan11@imperial.ac.uk, panweihit@gmail.com)
%
% ============= Version =============
%   1.0 (Sep ?, 2012)
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
    % The cost function where we minimise the weights
    % Right is norm term
    minimize    (lambda*norm( U(:,iter).*W, 1 )+ 0.5*sum((Dic* W-Output).^2) )
    %                 subject to
    %                           W.^2-ones(101,1)<=0;
    cvx_end
    
    w_estimate(:,iter)=W;
    WWW(:,iter)=W;
    Gamma(:,iter)=U(:,iter).^-1.*abs(W);
    Dic0=lambda*eye(M)+Dic*diag(Gamma(:,iter))*Dic';
    UU(:, iter)=diag(Dic'*(Dic0\Dic));
    U(:,iter+1)=abs(sqrt(UU(:, iter)));
    
    w_unpruned = w_estimate(:,iter);

    for i=1:N
        if   w_estimate(i,iter).^2/norm(w_estimate(:,iter))^2<delta
            w_estimate(i,iter)=0;
        end
    end
    
    %Termination of loop on convergence only if iter > 5
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

% Output only the last estimate, because amount of iterations is not known
w_out = w_estimate(:, iter);

% The final cost where the algorithm quits
cost = lambda*norm( U(:,iter).*w_out, 1 ) + 0.5*sum((Dic* w_out-Output).^2);
penalty = lambda*norm( U(:,iter).*W, 1 );
ols = 0.5*sum((Dic* w_out-Output).^2);
