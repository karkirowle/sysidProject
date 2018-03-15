% Returns the parameters originally set in GRN_dis
% Bence Halpern
function parameters = getGRNParameters

parameters.ga1 = 0.3;
parameters.ga2 = 0.4;
parameters.ga3 = 0.5;
parameters.ga4 = 0.2;
parameters.ga5 = 0.4;
parameters.ga6 = 0.6;

parameters.be1 = 1.4;
parameters.be2 = 1.5;
parameters.be3 = 1.6;

parameters.al1 = 4;
parameters.al2 = 3;
parameters.al3 = 5;

parameters.h = 4;


ga_diag = diag( [
    -parameters.ga1, ...
    -parameters.ga2, ...
    -parameters.ga3, ...
    -parameters.ga4, ...
    -parameters.ga5, ...
    -parameters.ga6, ...
    ]);

be_diag = diag([parameters.be1, parameters.be2, parameters.be3]);

% The coefficient matrix for the library functions
parameters.w_tru = [ ...,
    ga_diag + [zeros(3,3) be_diag;zeros(3,6)]; ...,
    [zeros(15,3);0 parameters.al2 0;0 0 parameters.al3;parameters.al1 0 0;zeros(30,3);] zeros(48,3)];

end