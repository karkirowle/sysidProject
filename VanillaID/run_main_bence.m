
clear all;
close all;



% Runs the repressilator simulation
repressilator

% The construction of the dictionary indicates the amount of assumed states
lambda = 0.01; 
MAXITER = 5;
states = 6;
time = 100;

Y = zeros(time-1,states);
for k=2:(time)
    
    
    Phi(k-1,:)=[ x_2(k-1,:), ...
               hill(x_2(k-1,:),1), ...
               hill(x_2(k-1,:),2), ... 
               hill(x_2(k-1,:),3), ... 
               hill(x_2(k-1,:),4)];
           
    
    Y(k-1,:)=  (x_2(k,:)-x_2(k-1,:)); 
end

%W_ours will be dimensions number of function * times output
% One nuance that I missed - the things that needs to be inputted is the
% difference, because it is a differential equation
w_ours =  tac_reconstruction(Y(:,1), Phi, lambda,MAXITER);

disp(w_ours);
%compare = [w_true, w_ours(:,end), w_lasso, w_ls]
