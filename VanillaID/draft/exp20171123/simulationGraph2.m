classdef simulationGraph2
    properties (Access = public)
      basisFunctions
      geneIdentifier
      weightList
      geneList
      numberOfGenes 
    end
    methods
        function obj = simulationGraph2(numberOfGenes)
            % each gene equation should have a packet of basis functions
            % These should store a packet of function
            obj.basisFunctions = {};
            obj.geneIdentifier = [];
            obj.weightList = [];
            obj.geneList = [];
            obj.numberOfGenes = numberOfGenes;
        end
        function obj = addBasisFunction(obj, fun, geneEquation, gene, weight)
            obj.basisFunctions(end+1) = {fun};
            obj.geneIdentifier(end+1) = geneEquation;
            obj.weightList(end+1) = weight;
            obj.geneList(end+1) = gene;
        end
       function derivativeSeries = runSimulation(obj, initialConditions, timeStep, noise, h)
           assert(length(initialConditions) == obj.numberOfGenes, ...
               'You havent provided the right amount initial conditions for the simulation!');
           derivativeSeries = zeros(timeStep-1, obj.numberOfGenes);
           results = zeros(timeStep,obj.numberOfGenes);
           results(1,:) = initialConditions;
        for i=2:timeStep 
           results(i,:) = results(i-1,:);
           for k = 1:obj.numberOfGenes
               geneIndex = find(obj.geneIdentifier == k);
               currentBases = obj.basisFunctions(geneIndex);
               for j = 1:length(currentBases)
                   currentIndex = geneIndex(j);
                   currentFunction = obj.basisFunctions(currentIndex);
                   currentWeight = obj.weightList(currentIndex);
                   currentData = results(i-1,obj.geneList(currentIndex));
                   results(i, k) = results(i, k) + h*currentWeight*currentFunction{1}(currentData);
               end
               results(i,k) = results(i,k) + noise*randn;
           end
            derivativeSeries(i-1,:)=(results(i,:)-results(i-1,:))/h;
        end
         end
        
    end
end
