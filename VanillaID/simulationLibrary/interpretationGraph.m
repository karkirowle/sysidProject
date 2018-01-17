classdef interpretationGraph
    properties (Access = public)
      basisFunctions
      geneIdentifier
      weightList
      geneList
      numberOfGenes 
    end
    methods
        function obj = interpretationGraph(numberOfGenes)
            obj.numberOfGenes = numberOfGenes;
            obj.basisFunctions = {};
            obj.weightList = [];
    
        end
        function obj = addBasisFunction(obj, fun)
            obj.basisFunctions(end+1) = {fun};
        end
        function Phi = constructDictionary(obj, timeSeries)
            timeSpan = size(timeSeries,1);
            states = size(timeSeries,2);
            numFunctions = length(obj.basisFunctions);
            for k=1:timeSpan
                for j=1:numFunctions
                    coeff1 = 1 + (j-1)*states;
                    coeff2 = j*states;
                    % This evaluates the function for all states/genes
                    Phi(k,coeff1:coeff2) = obj.basisFunctions{j}(timeSeries(k,:));
                end
            end
        end
        function [obj, w_estimate, cost] = reconstruct(obj, timeSeries, derivativeSeries)
            Phi = constructDictionary(obj,timeSeries);
            disp(size(Phi));
            disp(size(derivativeSeries));
            [estimate_temp, cost] = tac_reconstruction(derivativeSeries, Phi, 0, 5);
            w_estimate = estimate_temp(:,5);
        end
        function motifSeries = motifCalculation(obj, motifInterval, w_estimate, byWhich)
            counter = 0;
            numFunctions = length(obj.basisFunctions);
            coefficients = byWhich:obj.numberOfGenes:(obj.numberOfGenes*numFunctions);
            disp(coefficients);
                for j=1:numFunctions
                    counter = counter + 1;
                    motifSeries(counter,:) = w_estimate(coefficients(j))*obj.basisFunctions{j}(motifInterval);
                end
                motifSeries = sum(motifSeries,1);
        end
        
    end
end
