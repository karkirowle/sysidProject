% Euler method repressilator
function  [w_ours, error] = eulerRepressilator(state)
    samplingRate = 1;
    timePoints = 100;
    noise = 0.01;
    lambda = 0.01;
    MAXITER = 5;

    realState = 6;
    if (state > 6)
            X = zeros(timePoints,state);
        X(1,:) = [1,1,1,20,10,10, zeros(1,state-realState)];
    else
         X = zeros(timePoints,realState);
        X(1,:) = [1,1,1,20,10,10];
   
    end
    
    w_tru = zeros(54,realState);

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
               sHill(X(k-1,:),1,realState), ...
               sHill(X(k-1,:),2,realState), ... 
               sHill(X(k-1,:),3,realState), ... 
               sHill(X(k-1,:),4,realState)];
        Phi2(k-1,:) = [ X(k-1,:), ...
               sHill(X(k-1,:),1,state), ...
               sHill(X(k-1,:),2,state), ... 
               sHill(X(k-1,:),3,state), ... 
               sHill(X(k-1,:),4,state)];

        % Generate next time step of X 
        
        %disp(size(Phi(k-1,1:6)))
        %disp(size(w_tru))
        X(k,1:6)=X(k-1,1:6) + ...
            samplingRate*Phi(k-1,1:54)*w_tru+ ... % coeff determine which affect X
            + noise*randn(1,6); % added noise
        if (state > 6)
            X(k,7:state) = X(k-1,7:state) + randn;
        end
        % Take the derivative using a difference eq. and known sampling rate
        Y(k-1,:)=  (X(k,:)-X(k-1,:))/samplingRate; 
    end

    plot(1:100,X)
    xlabel("time")
    ylabel("concentration")
    title("Repressilator")

    w_ours =  tac_reconstruction(Y(:,1), Phi2, lambda,MAXITER);
    
    % Difference in the ratio of number of zero parameters to all
    % parameters (original should be higher and be approximated)
    error = length(find(w_ours(:,5)))/length(w_ours(:,5)) - ...,
        length(find(w_tru(:,1)))/length(w_tru(:,1));
end

