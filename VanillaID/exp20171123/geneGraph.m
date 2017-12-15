classdef geneGraph
    properties (Access = public)
       numberOfGenes
       odeList
    end
    methods
        function obj=geneGraph(amountOfGenes)
            obj.numberOfGenes = amountOfGenes;
        end
        function obj=degradation(obj, whichGene, degradationRate)
            obj.odeList{whichGene} =@(x)(obj.odeList{whichGene}(x) - degradationRate*x);
        end
        function obj=repression(obj, repressedBy, repressed, order, rate)
            obj.odeList{repressed} =@(x)(obj.odeList{whichGene}(x) - ... 
                (1./(1+x.^h)));
        end
        function obj=activation(obj, activatedBy, activated, order, rate)
            obj.odeList{whichGene} =@(x)(obj.odeList{whichGene}(x) -  ...
                x.^h/(1+x.^h));
        end
        function results=runSimulation(simulationTime, timeStep, initialConditions)
            results = zeros(obj.numberOfGenes, simulationTime/timeStep);
            results(:,1) = initialConditions;
            for i=2:(simulationTime/timeStep)
                for j=1:obj.numberOfGenes
                    results(j,i) = results(j,i-1) + ...
                        timeStep*obj.odeList{j}(results(j,i-1));
                end
            end
        end
    end
   
end