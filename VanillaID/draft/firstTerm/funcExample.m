% example Function Handle approach

func1 = @(x) 1 ./ (1 + x);
func2 = @(x) 1 ./  (1 + x^2);

totalFunc = {func1; func2};
lvar = length(totalFunc);
totalFunc = {totalFunc{lvar,1}; func1};
disp(totalFunc);
