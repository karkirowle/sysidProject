% Repressilator consistency experiment wrapper

% Experimental setting 1
% We increase the number of observable states to the original states and
% see if the error decreases, thus proving that the system can correctly
% identify how many measurements are needed?

stateToTest = 15;
error = zeros(1,stateToTest);
for state=1:stateToTest
    [wthrow, error(state)] = eulerRepressilator(state);
end

plot(1:stateToTest,error);
xlabel("number of states used");
ylabel("difference of sparsity ratios");
title("Error curve");