% Fisher Information for basis function

x = 1:10;

basis1 = 1./(1+x);
basis2 = 1./(1+x).^2;

% One data point 2 *1 1*2

for i=1:2
    x = [basis1(i); basis2(i)];
    fisher = x * x';
    fisherDet(i) = det(fisher);
end

x = [basis1(1:2); basis2(1:2)];
fisher = x * x';
fisherDet2 = det(fisher);