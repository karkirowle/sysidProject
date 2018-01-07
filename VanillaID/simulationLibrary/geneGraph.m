classdef geneGraph
    properties (Access = public)
        simulationGraph
    end
    methods
        function obj = geneGraph(numberOfGenes)
            obj.simulationGraph = simGraph(numberOfGenes);
        end
        function obj = degradation(obj, whichGene, degradationRate)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) x, whichGene, whichGene, degradationRate);
        end
        function obj = activation(obj, activated, activator, activationRate, hillOrder)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) x./(1+x).^hillOrder, activated, activator, activationRate);
        end
        function obj = repression(obj, repressed, repressor, repressionRate, hillOrder)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) 1./(1+x).^hillOrder, repressed, repressor, repressionRate);
        end
        function [derivativeSeries, timeSeries] = runSimulation(obj, initialConditions, timeStep, noise, h)
            [derivativeSeries, timeSeries] = obj.simulationGraph.runSimulation(initialConditions, timeStep, noise, h);
        end
    end
end

    
    