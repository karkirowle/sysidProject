

function hillFunction = sHill( x,h,s )

% This Hill function subroutine does the following thing:
% Accepts a state vector with a length of s
% Accepts h, the order of the hill function desired


% Evaluates the Hill function for all states bot in activating and in
% repressing form

% x - input state vector with length of 6
% h - order of Hill Function

% Preallocation
hillFunction= zeros(1,2*s);

% Generate function
for i=1:2*s
    if i<=s
        hillFunction(i)=1/(1+x(i)^h);
    else
        hillFunction(i)=x(i-s)^h/(1+x(i-s)^h);
    end
end
end

