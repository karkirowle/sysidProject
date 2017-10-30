% Repressilator consistency experiment wrapper

% Experimental setting 1
% We increase the number of observable states to the original states and
% see if the error decreases, thus proving that the system can correctly
% identify how many measurements are needed?

error = zeros(1,6);
for state=1:6
    [wthrow, error(state)] = eulerRepressilator(state);
end

plot(1:6,error);
xlabel("number of states used");
ylabel("difference of sparsity ratios");
title("Error curve");