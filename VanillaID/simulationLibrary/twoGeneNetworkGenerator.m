function [timeSeries, derivativeSeries, coeff] =twoGeneNetworkGenerator(numNet, time, noise, h)
% numNets - number of networks to generate
% seed - option to use a constant seed (reproducibility aims)

timeSeries = zeros(numNet, time,2);
derivativeSeries = zeros(numNet, time-1, 2);
coeff = zeros(4,numNet);

for i=1:numNet
    sim = geneGraph(2);
    
    coeff(1,i) = -0.9*rand;
    coeff(2,i) = -0.9*rand;
    coeff(3,i) = 100*rand;
    coeff(4,i) = 100*rand;
    sim = sim.degradation(1, coeff(1,i));
    sim = sim.degradation(2, coeff(2,i));
    sim = sim.repression(1,2, coeff(3,i),4);
    sim = sim.repression(2,1, coeff(4,i) ,4);
    
    initialConditions = 5 + 100*rand(1,2);
    
    [derivativeSeriesTemp, timeSeriesTemp] = ...,
        sim.runSimulation(initialConditions, time, noise, h);

    timeSeries(i,:,:) = timeSeriesTemp;
    derivativeSeries(i,:,:) = derivativeSeriesTemp;
    
end

end
