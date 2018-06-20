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
        function Phi = constructDictionary(obj, timeSeries, bias)
            timeSpan = size(timeSeries,1);
            states = size(timeSeries,2); % TODO: This should be equal to number of genes
            numFunctions = length(obj.basisFunctions);
            for k=1:timeSpan
                for j=1:numFunctions
                    coeff1 = 1 + (j-1)*states;
                    coeff2 = j*states;
                    %This evaluates the function for all states/genes
                    Phi(k,coeff1:coeff2) = obj.basisFunctions{j}(timeSeries(k,:));
                end
                if (bias)
                    Phi(k, numFunctions*states + 1) = 1;
                end
            end
        end
        function [obj, w_estimate, cost] = reconstruct(obj, timeSeries, derivativeSeries, lambda)
            Phi = constructDictionary(obj,timeSeries, false);
            [estimate_temp, cost, ~, ~, ~, ~, ~] = tac_reconstruction(derivativeSeries, Phi, lambda, 5);
            w_estimate = estimate_temp(:,5);
        end
        function [obj, w_estimate, cost, w_unprunedOut, penalty, ols] = reconstructUnpruned(obj, timeSeries, derivativeSeries, lambda, bias)
            [w_estimate, cost, w_unprunedOut, penalty, ols, ~, ~]  = ...,
                reconstructSetIter(timeSeries, derivativeSeries, ...,
                lambda, bias, 5);
        end
        function [obj, w_estimate, cost, w_unprunedOut, penalty, ols, convergenceGamma] = reconstructSetIter(obj, timeSeries, derivativeSeries, lambda, bias, iter)
            Phi = constructDictionary(obj,timeSeries, bias);
            [w_estimate, cost, w_unpruned, penalty, ols, convergenceGamma, ~] = ...,
                tac_reconstruction(derivativeSeries, Phi, lambda, iter);
            w_unprunedOut = w_unpruned;
        end
        function [obj, w_estimate, cost, w_unprunedOut, penalty, ols, convergenceGamma, Gamma] = reconstructGamma(obj, timeSeries, derivativeSeries, lambda, bias, iter)
            Phi = constructDictionary(obj,timeSeries, bias);
            [w_estimate, cost, w_unpruned, penalty, ols, convergenceGamma, Gamma] = ...,
                tac_reconstruction(derivativeSeries, Phi, lambda, iter);
            w_unprunedOut = w_unpruned;
        end
        function motifSeries = motifCalculation(obj, motifInterval, w_estimate, byWhich)
            counter = 0;
            numFunctions = length(obj.basisFunctions);
            coefficients = byWhich:obj.numberOfGenes:(obj.numberOfGenes*numFunctions);
            for j=2:numFunctions
                counter = counter + 1;
                motifSeries(counter,:) = w_estimate(coefficients(j))*obj.basisFunctions{j}(motifInterval);
            end
            motifSeries = sum(motifSeries,1);
        end
        
    end
end
