

function hillFunction = hill( x,h )

% This Hill function subroutine does the following thing:
% Accepts a state vector with a length of 6
% Accepts h, the order of the hill function desired


% Evaluates the Hill function for all states bot in activating and in
% repressing form

% x - input state vector with length of 6
% h - order of Hill Function

% Preallocation
hillFunction= zeros(1,12);

% Generate function
for i=1:12
    if i<=6
        hillFunction(i)=1/(1+x(i)^h);
    else
        hillFunction(i)=x(i-6)^h/(1+x(i-6)^h);
    end
end
end

