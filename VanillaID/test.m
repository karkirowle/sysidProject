% Sample code - Bence Halpern

% GRN dis
whichState = 1; % You can select the returned state varaiable from 1-6
noise = 0.01; % Noise applied to the Gaussian
parameters = getGRNParameters;
[a, b, c] = GRN_dis(whichState, noise, parameters);