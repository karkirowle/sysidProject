% Three gene net for reconstruction method
function  [w_ours, error] = threegenenet(state, which)
% state - the number of states to use for predicting (a number greater
% equal than 1)
% which - which state to predict (a number from 1-6)

samplingRate = 1;
timePoints = 100;
noise = 0.01;
lambda = 0.1;
MAXITER = 5;
realState = 3;
h = 4;

if (state > 3)
    X = zeros(timePoints,state);
    X(1,:) = [20,10,10, zeros(1,state-realState)];
else
    X = zeros(timePoints,realState);
    X(1,:) = [20,10,10];
    
end

functionNumber = realState + h*2*realState;
w_tru = zeros(realState + h*2*realState,realState);
% 25,26,27 repress g1-g2-g3
% 22,23,24 activate g1-g2-g3
% Gene 1 - inhibited by gene 3, degradation
w_tru(1,1) = -0.20;
w_tru(27,1) = 0.3;
% Gene 2 - inhibited by gene 1, degradation
w_tru(2,2) = -0.20;
w_tru(25,2) = 0.3;
% Gene 3 - activated by gene 2, degradation
w_tru(3,3) = -0.20;
w_tru(23,3) = 0.3;


for k=2:timePoints
    Phi(k-1,:) = [ X(k-1,:), ...
        sHill(X(k-1,:),1,realState), ...
        sHill(X(k-1,:),2,realState), ...
        sHill(X(k-1,:),3,realState), ...
        sHill(X(k-1,:),4,realState)];
    % Here we need to modify that not for all generated, but rather just
    % the number observed
    Phi2(k-1,:) = [ X(k-1,1:state), ...
        sHill(X(k-1,1:state),1,state), ...
        sHill(X(k-1,1:state),2,state), ...
        sHill(X(k-1,1:state),3,state), ...
        sHill(X(k-1,1:state),4,state)];
    
    % Generate next time step of X
    
%     disp(size(Phi(k-1,1:functionNumber)))
%     disp(size(w_tru))
    X(k,1:realState)=X(k-1,1:realState) + ...
        samplingRate*Phi(k-1,1:functionNumber)*w_tru+ ... % coeff determine which affect X
        + noise*randn(1,3); % added noise
    if (state > realState)
        X(k,(realState+1):state) = X(k-1,(realState+1):state) + noise*randn(1,3);
    end
    % Take the derivative using a difference eq. and known sampling rate
    Y(k-1,:)=  (X(k,:)-X(k-1,:))/samplingRate;
end

figure(1)
plot(Y);
figure(2);
plot(1:100,X)
xlabel("time")
ylabel("concentration")
title("Repressilator")
legend('gene 1', 'gene 2', 'gene 3')
w_ours =  tac_reconstruction(Y(:,which), Phi2, lambda,MAXITER);

% Difference in the ratio of number of zero parameters to all
% parameters (original should be higher and be approximated)
error = length(find(w_ours(:,5)))/length(w_ours(:,5)) - ...,
    length(find(w_tru(:,which)))/length(w_tru(:,which));
end

