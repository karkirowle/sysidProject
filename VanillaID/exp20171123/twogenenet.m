% Two gene net for reconstruction method
% Aim: Gene 1 represses Gene 
function  [w_ours, w_tru, w_compare, funcListEnd] = twogenenet(state, which)
% state - the number of measurement sets to to use for predicting 
% (a number greater  equal than 1)
% which - which state to predict (a number from 1-6)

% --------------------------- Parameters ----------------------------------

samplingRate = 1;
timePoints = 100;
noise = 0;
lambda = 0;
MAXITER = 5;
realState = 2;
h = 4;
initialConditions = [20,10]; 

% Interpretation
timeLag = 0;


% ------------------------- Initial conditions ----------------------------

% If more measurement states are supposed generate extra initial conditions

if (state > realState)
    X = zeros(timePoints,state);
    % TODO: Maybe ones would be a better choice?
    X(1:(1+timeLag),:) = repmat( [initialConditions, ...
            zeros(1,state-realState)], timeLag+1, 1);
else
    %X = zeros(timePoints,state);
    X(1:(1+timeLag),:) = repmat(initialConditions, timeLag+1, 1);
end

disp("Start of X")
disp(X)
disp("End of initial X");

% ------------------------- Simulation graph ------------------------------

% Simulation settings
simulation = simulationGraph(2,h,0,1,false);
simulation = simulation.degradation(1,0.2);
simulation = simulation.degradation(2,0.2);
simulation = simulation.repression(1,2, 0.5, h);
simulation = simulation.repression(2,1, 0.5, h);
w_tru = simulation.weightMatrix;
functionNumber = size(w_tru,1);

% Interpretation graph -> expecting state many nodes with order of 4
%interpretation = simulationGraph(state,4,timeLag,1,true);
% coincides with simulation when timeLag = 0, bias = false, state = 1
interpretation = simulationGraph(1,h,timeLag,1,false); % here i modi

% ------------------------- Euler simulation ------------------------------

for k=(2+timeLag):timePoints
    % Simulation dictionary functions
    [Phi(k-1,:), funcList1] = simulation.getNextDictionary(X(k-1,:));
    % feed 20 10 10 in
    % Identification dictionary functions
    % Uses only the number observed
    % feed 20 20 10 two times in, but one gene so that 20;20
    [Phi2(k-1,:), funcList2] = interpretation.getNextDictionary(X((k-1-timeLag) ...
        :(k-1),1:state));
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
disp(size(Y(:, which)));
disp(size(Phi2));
[w_ours,cost] =  tac_reconstruction(Y(:,which), Phi2, lambda,MAXITER);


% ------------------------ Error analysis ---------------------------------
disp(['Loss value at end', num2str(cost)]);

% Compare reconstruction errors if equal (L1-norm) 
if (length(w_tru) == size(w_ours,1))
    errorL1 = sum(norm(w_ours(:,MAXITER)-w_tru,1))/length(w_tru);
    disp(['L1 norm error:', num2str(errorL1)]);
    w_compare(:,1) = w_ours(:,MAXITER);
    w_compare(:,2) = w_tru(:,which);
else
    w_compare = []; 
end


% ---------------------- Signal reconstruction ----------------------------

% Graphical comparison of derivative of measurements
figure(1)
subplot(1,2,1);
plot(2:timePoints,Y, 'LineWidth', 1.5); % TOCHECK: Starting point for this
title('Derivative measurement curves');
xlabel('time');
ylabel('concentration');
box off;

linearApprox= Phi2 *w_ours(:,MAXITER);
subplot(1,2,2);
plot(2:timePoints, linearApprox, 'LineWidth', 1.5, 'Color', 'red');
title('Selected genes approximation');
xlabel('time');
ylabel('concentration');
box off;

% ----------------------- Motif interpretation ----------------------------

% Interpretation graph

interpretation.weightMatrix = w_ours(:,MAXITER); 
figure(2)
funcListEnd = interpretation.plotMotif(0:0.001:5,which);
title('Reconstructed motif of two-gene toy example')
xlabel('Increasing concentration of gene 1')
ylabel('Concentration effect on itself');
box off;

funcListEnd=funcListEnd.'


end

