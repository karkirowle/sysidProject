% Euler method repressilator

samplingRate = 1;
timePoints = 100;
noise = 0.01;
state = 6;
X = zeros(timePoints,state);
X(1,:) = [1,1,1,20,10,10];

w_tru = zeros(54,state);

% Gene 1 - repressed by protein3, degradation
w_tru(1,1) = -0.20;
w_tru(48,1) = 10;
% Gene 2 - repressed by protein1, degradation
w_tru(2,2) = -0.20;
w_tru(46,2) = 10;
% Gene 3 - repressed by protein2, degradation
w_tru(3,3) = -0.20;
w_tru(47,3) = 10;
% Protein 1 - transcribed by gene, degradation
w_tru(4,4) = -0.20;
w_tru(1,4) = 0.20;
% Protein 2 - transcribed by gene, degradation
w_tru(5,5) = -0.20;
w_tru(2,5) = 0.20;
% Protein 1 - transcribed by gene, degradation
w_tru(6,6) = -0.20;
w_tru(3,6) = 0.20;

for k=2:timePoints 
    Phi(k-1,:) = [ X(k-1,:), ...
           sHill(X(k-1,:),1,state), ...
           sHill(X(k-1,:),2,state), ... 
           sHill(X(k-1,:),3,state), ... 
           sHill(X(k-1,:),4,state)];
    Phi2(k-1,:) = [ X(k-1,:), ...
           sHill(X(k-1,:),1,7), ...
           sHill(X(k-1,:),2,7), ... 
           sHill(X(k-1,:),3,7), ... 
           sHill(X(k-1,:),4,7)];
       
    % Generate next time step of X 
    
    X(k,:)=X(k-1,:) + ...
        samplingRate*Phi(k-1,:)*w_tru+ ... % coeff determine which affect X
        + noise*randn(1,6); % added noise
    
    % Take the derivative using a difference eq. and known sampling rate
    Y(k-1,:)=  (X(k,:)-X(k-1,:))/samplingRate; 
end

plot(1:100,X)

w_ours =  tac_reconstruction(Y(:,1), Phi2, lambda,MAXITER);
