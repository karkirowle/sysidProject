% Let's write a simple repressilator

clc;
clear;
close all;

a = 10;
b = 10;
c = 10;
%e = 0.22;
d = 0.20;

[t,x] = ode45(@(t,x)repressilatorFun(x,a,b,c,d),[0 100],[1,1,1,20,10,10]);

% I assume the output function needs to be interpolated for hundred points
% Interpolation step
queryPoints = 1:1:100;
x_2 = interp1(t,x,queryPoints);
plot(queryPoints,x_2)
title("Repressilator")
xlabel("time")
ylabel("concentration")


function dydt = repressilatorFun(x,a,b,c,d)
    noise = 0.01;
    gene1 = x(1);
    gene2 = x(2);
    gene3 = x(3);
    p1 = x(4);
    p2 = x(5);
    p3 = x(6);
    
    dydt = zeros(6,1);
    dydt(1) = a./(1+p3.^4) - d.*gene1 + noise*randn;
    dydt(2) = b./(1+p1.^4) - d.*gene2 + noise*randn;
    dydt(3) = c./(1+p2.^4) - d.*gene3 + noise*randn;
    dydt(4) = d.*(gene1 - p1) + noise*randn;
    dydt(5) = d.*(gene2 - p2) + noise*randn;
    dydt(6) = d.*(gene3 - p3) + noise*randn;
end