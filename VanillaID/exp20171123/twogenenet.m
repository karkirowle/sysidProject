% Two gene net for reconstruction method
% Aim: Gene 1 represses Gene 2

function  [w_ours] = twogenenet(state, which)
% state - the number of measurement sets to to use for predicting 
% (a number greater  equal than 1)
% which - which state to predict (a number from 1-6)

% --------------------------- Parameters ----------------------------------

samplingRate = 1;
timePoints = 100;
noise = 0.01;
lambda = 0.1;
MAXITER = 5;
realState = 2;
h = 4;
initialConditions = [20,10,10];

% ------------------------- Initial conditions ----------------------------

% If more measurement states are supposed generate extra initial conditions

if (state > realState)
    X = zeros(timePoints,state);
    % TODO: Maybe ones would be a better choice?
    X(1,:) = [initialConditions, zeros(1,state-realState)];
else
    X = initialConditions;
end

% ------------------------- Simulation graph ------------------------------

% Simulation settings
simulation = simulationGraph(2,4);
simulation = simulation.degradation(1,0.2);
simulation = simulation.degradation(2,0.2);
simulation = simulation.repression(1,2, 0.3, 4);
simulation = simulation.repression(2,1, 0.3, 4);
w_tru = simulation.weightMatrix;
functionNumber = size(w_tru,1);

% ------------------------- Euler simulation ------------------------------

for k=2:timePoints
    % Simulation dictionary functions
    Phi(k-1,:) = [ X(k-1,:), ...
        sHill(X(k-1,:),1,realState), ...
        sHill(X(k-1,:),2,realState), ...
        sHill(X(k-1,:),3,realState), ...
        sHill(X(k-1,:),4,realState)];
    
    % Identification dictionary functions
    % Uses only the number observed
    Phi2(k-1,:) = [ X(k-1,1:state), ...
        sHill(X(k-1,1:state),1,state), ...
        sHill(X(k-1,1:state),2,state), ...
        sHill(X(k-1,1:state),3,state), ...
        sHill(X(k-1,1:state),4,state)];
    
    % Generate next time step of X
    X(k,1:realState)=X(k-1,1:realState) + ...
        samplingRate*Phi(k-1,1:functionNumber)*w_tru + ... 
        + noise*randn(1,realState); 
    
    % Generate trajectories from noise when no measurement available for
    % additional state
    if (state > realState)
        X(k,(realState+1):state) = X(k-1,(realState+1):state) ...
            + noise*randn(1,(realState+1):state);
    end
    
    % Take the derivative using a difference eq. and known sampling rate
    Y(k-1,:)=  (X(k,:)-X(k-1,:))/samplingRate;
end

% ------------------------- Reconstruction --------------------------------

w_ours =  tac_reconstruction(Y(:,which), Phi2, lambda,MAXITER);
disp(w_ours);

end

