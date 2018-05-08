
% Genetic Regulatory Network distribution
function [y, A,w_true] = GRN_dis(whichState,noise,param)

% whichState - "which state variable you are interested in
% noise - Factor for multiplying Gaussian Noise with mean 0
% param - the parameters which are originally hardcoded in the model
% A - some kind of functions used ???

% state - 6 states are HARD-CODED
% h - 4 is the order, HARD-CODED
% w=rand(6+h*12,6);
% w_tru=[diag([-parameters.ga1,-parameters.ga2,-parameters.ga3,-parameters.ga4,-parameters.ga5,-parameters.ga6])+[zeros(3,3) diag([parameters.be1,parameters.be2,parameters.be3]);zeros(3,6)];[zeros(45,3);0 parameters.al2 0;0 0 parameters.al3;parameters.al1 0 0] zeros(48,3)];

% Getting ground truth matrix from parameters
w_tru = param.w_tru;

% Number of time points
timePoints = 101;

% Varables acronym for the derivate vector size
M = timePoints - 1;

% Records of states: T rows, 6 states, T*6
X = zeros(timePoints,6);

% Set initial states from a normal distribution btwn 0-1 for all 6
X(1,:)=rand(1,6);

% Preallocation
% 6 + 4*12 = 54 -> numbers of states + orders of Hill used * 2 * number of
% states both for activation and repression

Phi=zeros(M,6+param.h*12);
Y=zeros(M,6);

samplingRate = 1;

% Iterate through each time step
for k=2:timePoints 
    %real time from 0:1:50
    %call hill function:H
    
    % Phi contains substituted results for all states into activated and
    % repressed Hill functions 
    Phi(k-1,:)=[ X(k-1,:), ...
           hill(X(k-1,:),1), ...
           hill(X(k-1,:),2), ... 
           hill(X(k-1,:),3), ... 
           hill(X(k-1,:),4)];
    
    % Generate next time step of X 
    X(k,:)=X(k-1,:) + ...
        samplingRate*Phi(k-1,:)*w_tru+ ... % coeff determine which affect X
        + noise*randn(1,6); % added noise
    
    % Take the derivative using a difference eq. and known sampling rate
    Y(k-1,:)=  (X(k,:)-X(k-1,:))/samplingRate; 
end


% Output the state variable the user is interestedin
y = Y(:,whichState);

% Rename phi to A
A = Phi;

% Acquire the desired column from the ground thruth matrix 
w_true = w_tru(:,whichState);

% Calculate and displays the errors which are due to the added noise
% I assume the purpose of this is to calculate SNR maybe...?
err = norm(y-A*w_true);
disp(err);

% Creates a plot from X
figure; 
plot(X); 
