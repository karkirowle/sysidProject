classdef geneGraph
    properties (Access = public)
        simulationGraph
        digraph
        highlight
    end
    methods
        function obj = geneGraph(numberOfGenes)
            obj.simulationGraph = simGraph(numberOfGenes);
            obj.digraph = digraph();
            for i=1:numberOfGenes
                obj.digraph = addnode(obj.digraph,{['Gene ' num2str(i)]});
            end
            obj.highlight = 0;
        end
        function obj = degradation(obj, whichGene, degradationRate)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) x, whichGene, whichGene, degradationRate);
        end
        function obj = activation(obj, activated, activator, activationRate, hillOrder)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) x./(1+x).^hillOrder, activated, activator, activationRate);
            obj.digraph = addedge(obj.digraph,{['Gene ' num2str(activator)]}, {['Gene ' num2str(activated)]}, activationRate);
            obj.highlight = [obj.highlight 0];
        end
        function obj = repression(obj, repressed, repressor, repressionRate, hillOrder)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) 1./(1+x).^hillOrder, repressed, repressor, repressionRate);
            obj.digraph = addedge(obj.digraph,{['Gene ' num2str(repressor)]}, {['Gene ' num2str(repressed)]}, repressionRate);
            obj.highlight = [obj.highlight 1];
        end
        function [derivativeSeries, timeSeries] = runSimulation(obj, initialConditions, timeStep, noise, h)
            [derivativeSeries, timeSeries] = obj.simulationGraph.runSimulation(initialConditions, timeStep, noise, h);
        end
        function [derivativeSeries, results, timePoints] = runRungeKutta(obj, initialConditions, noise, tspan)
            [derivativeSeries, results, timePoints] = obj.simulationGraph.runRungeKutta(initialConditions, noise, tspan);
        end
        function showGraph(obj, ~)
            h = plot(obj.digraph, 'EdgeLabel', obj.digraph.Edges.Weight);
            for i=2:length(obj.highlight)
                if (obj.highlight(i) == 1)
                    geneNodes= table2array(obj.digraph.Edges(i-1,1));
                    highlight(h,{geneNodes{1}},{geneNodes{2}},'EdgeColor','r');
                end
            end
        end
    end
end

    
    