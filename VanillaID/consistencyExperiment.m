% Repressilator consistency experiment wrapper

% Experimental setting 1
% We increase the number of observable states to the original states and
% see if the error decreases, thus proving that the system can correctly
% identify how many measurements are needed?

% stateToTest = 15;
% allStates = 6;
% error = zeros(allStates,stateToTest);
% for measurement=1:allStates
%     for state=1:stateToTest
%         [wthrow, error(measurement,state)] = eulerRepressilator(state, ...
%             measurement);
%     end
% end
% 
% 
% 
% plot(1:stateToTest,error(1,:));
% hold on;
% plot(1:stateToTest,error(2,:));
% hold on;
% plot(1:stateToTest,error(3,:));
% hold on;
% plot(1:stateToTest,error(4,:));
% hold on;
% plot(1:stateToTest,error(5,:));
% hold on;
% plot(1:stateToTest,error(6,:));
% legend('state 1','state2','state3','state4','state5','state6')
% xlabel("number of states used");
% ylabel("difference of sparsity ratios");
% title("Error curve");