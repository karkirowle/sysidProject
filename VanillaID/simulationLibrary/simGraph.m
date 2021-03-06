classdef simGraph
    properties (Access = public)
      basisFunctions
      geneIdentifier
      weightList
      geneList
      numberOfGenes 
    end
    methods
        function obj = simGraph(numberOfGenes)
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
        function [derivativeSeries, results, timePoints] = runRungeKutta(obj, initialConditions, noise, tspan)
 
            assert(length(initialConditions) == obj.numberOfGenes, ...
                'You havent provided the right amount initial conditions for the simulation!');
            [t,x]=ode45(@(t,x) rungeKuttaHelper(obj,t,x,noise),tspan,initialConditions,noise);
            results = x;
            timePoints = t;
            derivativeSeries = derivativeEvaluator(obj,results);
        end
        function [f] = rungeKuttaHelper(obj,~,x,noise)
            f = zeros(obj.numberOfGenes,1);
            for k = 1:obj.numberOfGenes
                geneIndex = find(obj.geneIdentifier == k);
                currentBases = obj.basisFunctions(geneIndex);
                for j = 1:length(currentBases)
                    currentIndex = geneIndex(j);
                    currentFunction = obj.basisFunctions(currentIndex);
                    currentWeight = obj.weightList(currentIndex);
                    currentData = x(obj.geneList(currentIndex));
                    f(k) = f(k) + currentWeight*currentFunction{1}(currentData);
                end
                f(k) = f(k) + noise*randn;
            end
        end
        function derivativeSeries = derivativeEvaluator(obj,timeSeries)
            derivativeSeries = zeros(size(timeSeries));
           for k = 1:obj.numberOfGenes
               geneIndex = find(obj.geneIdentifier == k);
               currentBases = obj.basisFunctions(geneIndex);
               for j = 1:length(currentBases)
                   currentIndex = geneIndex(j);
                   currentFunction = obj.basisFunctions(currentIndex);
                   currentWeight = obj.weightList(currentIndex);
                   currentData = timeSeries(:,obj.geneList(currentIndex));
                   derivativeSeries(:, k) = derivativeSeries(:,k) + currentWeight*currentFunction{1}(currentData);
               end
           end
        end

       function [derivativeSeries, results] = runSimulation(obj, initialConditions, timeStep, noise, h)
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
