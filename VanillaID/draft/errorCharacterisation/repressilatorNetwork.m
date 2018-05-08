% Housekeeping
clc;
clear;
close all;


% Creating eight gene repressing topology -> control example

% Sampling parameters from an uniform distribution interval [0,1]
% 8x1 = 8*1 1*4

% Creating ten realisations 
figure;
for i=1:10
parameters = (0.8 + 0.4*rand(8,1)) * [40, 1, 3, 0.5, 1];
initialConditions = rand(1,8);
% Repressilator
[t,x] = ode45(@(t,x) repressilatorHelper(x,parameters) , [0 1000], ...,
    initialConditions);
hold on;
plot(x(:,4));

end


function dxdt = repressilatorHelper(x, p)
% Parameters is the initial conditions as defined in the paper

dxdt = zeros(8,1);

dxdt(1) = p(1,1)/(1 + x(8)^1) - p(1,5)*x(1);
for i=2:8
   dxdt(i) = p(i,1)/(1+x(i-1)^1) -  p(i,5)*x(i); 
end

end
% 
% function dxdt = repressilatorHelper(x, p)
% % Parameters is the initial conditions as defined in the paper
% 
% dxdt = zeros(8,1);
% 
% dxdt(1) = p(1,1)/(p(1,2)^p(1,3) + x(8)^p(1,3)) + p(1,4) - p(1,5)*x(1);
% for i=2:8
%    dxdt(i) = p(i,1)/(p(i,2)^p(i,3)+x(i-1)^p(i,3)) + p(i,4) -  p(i,5)*x(i); 
% end
% 
% end

