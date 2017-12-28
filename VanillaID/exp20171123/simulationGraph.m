classdef simulationGraph
   properties (Access = public)
       numberOfGenes
       hillOrder
       weightMatrix
       timeLag % 0 is no timeLag
       amountOfHill
       amountOfNonHill
       polynomialOrder
       hillPerGene
       dictionaryLength
       functionPerHill
       bias
   end
   methods
       function obj = simulationGraph(numberOfGenes, ...
               hillOrder, timeLag, polynomialOrder, bias)
           % NOTE: timeLag is only for interpretations
           obj.numberOfGenes = numberOfGenes;
           obj.hillOrder = hillOrder;
           obj.timeLag = timeLag; 
           obj.polynomialOrder = polynomialOrder;
           obj.amountOfNonHill = numberOfGenes*(timeLag+1)*polynomialOrder;
           obj.bias = bias;
           if (bias)
            obj.amountOfNonHill = obj.amountOfNonHill + 1;
           end
           obj.amountOfHill = (timeLag+1)*hillOrder*2*numberOfGenes;
           obj.hillPerGene = obj.amountOfHill/obj.numberOfGenes;
           obj.dictionaryLength = obj.amountOfNonHill + obj.amountOfHill;
           obj.functionPerHill = (obj.dictionaryLength ...
               - obj.amountOfNonHill)/obj.hillOrder;
           obj.weightMatrix = zeros(obj.dictionaryLength,obj.numberOfGenes);
       end
       function obj = degradation(obj, degraded, amount)
           obj.weightMatrix(degraded,degraded) = -amount;
       end
       function obj = repression(obj, repressedBy, repressed, amount, order)
           % Total linear states + total previous order states + repressor
           % gene numbers
           row = obj.amountOfNonHill ...
               + (order-1)*obj.numberOfGenes*2 ...
               + repressedBy;
           
           obj.weightMatrix(row,repressed) = amount;
       end
       function obj = activation(obj, activatedBy, activated, amount, order)
           % Total linear states + total previous order states + repressor
           % states + total activation numbers
           row = obj.amountOfNonHill ...
               + (order-1)*obj.numberOfGenes*2 ...
               + obj.numberOfGenes ...
               + activatedBy;
           obj.weightMatrix(row,activated) = amount;
       end
       function [Phi, funcList] = getNextDictionary(obj, X)
           % Construct next dictionary step from given X at lag/cur t
           % Interpretation boolean flag is for the bias term
           counter = 0;

           % Add bias term to dictionary if this is an intepretation graph
           if (obj.bias)
               
               Phi = 1;
               counter = counter+1;
               funcList{counter} = 'bias';
           else
               
               Phi= [];
          
           end
           
           % Add polynomial terms to dictionary
           for j =1:(obj.timeLag+1) % TODO: Consider again what this timelag does, consider order of fors compared to Hill
               for k=1:obj.polynomialOrder
                   counter = counter+1;
                   funcList{counter} = ['polynom', num2str(k)];
                   Phi = [Phi, X(j,:).^k]; % Linear X(1,:)
               end
           end
           
         
           % Add Hill terms to dictionary
           % Ordered by delay descending and order ascending
           % So H1 H(t-1)1 H2 H2(t-1)
           for i=1:obj.hillOrder
               for j=1:(obj.timeLag+1)
                   counter = counter+1;
                   funcList{counter} = ['hill' num2str(i)];
                   Phi = [Phi, sHill(X(j,:),i,obj.numberOfGenes)];
               end
           end
       end
       function funcList = plotMotif(obj, interval, gene)
           % Uses the weightmatrix and a given interval to plot motif
           counter = 0;
           % Adding bias term
           if (obj.bias)
               counter = counter + 1;
               funcList{counter} = 'bias';
            hillTerm(counter,:) = obj.weightMatrix(1,gene) .* interval;
            
           end
           % Adding polynomial terms
           for l=1:(obj.timeLag+1)
               for k=1:(obj.polynomialOrder)
               % Interval selection
               currentInterval = zeros(size(interval));
               currentInterval(l:end) = interval(l:end);
               % Coefficient determination
               coeff = obj.bias + gene + obj.numberOfGenes*(k-1) + ...
                   obj.numberOfGenes*(l-1);
               disp(coeff)
               counter = counter + 1;
               funcList{counter} = ['polynomial', ...
                   'order:' num2str(k), ...
                   'timeLag:' num2str(l), ...
                   'coeff:' num2str(coeff)];
               hillTerm(counter,:) = obj.weightMatrix(coeff,1) .* ...
                   currentInterval.^k;
               end
           end
           % Adding Hill terms
           for h=1:obj.hillOrder
               for l=1:(obj.timeLag+1)
                   currentInterval = zeros(size(interval));
                   currentInterval(l:end) = interval(l:end);
                   for k=1:2
                       if (mod(k,2) == 1)
                           coeff = obj.bias + ...
                               obj.amountOfNonHill + ...
                               (h-1)*obj.functionPerHill + ...
                               (l-1)*(obj.functionPerHill/(obj.timeLag+1)) + ...
                               gene;
                           counter = counter + 1;
                           funcList{counter} = ['hill repressive ', ...
                               'order: ' num2str(k), ...
                               'timeLag: ' num2str(l), ...
                               'coeff: ', num2str(coeff)];
                            hillTerm(counter,:) = obj.weightMatrix(coeff ...
                                ,1) .* ...
                                (1./(1+interval.^h));
                       end
                       if (mod(k,2) == 0) % else is also valid
                           coeff = obj.bias + ...
                               obj.amountOfNonHill + ...
                               (h-1)*obj.functionPerHill + ...
                               (l-1)*(obj.functionPerHill/(obj.timeLag+1)) + ...
                               obj.numberOfGenes + ...
                               gene;
                        
                           counter = counter + 1;
                      funcList{counter} = ['hill activation', ...
                          'order:' num2str(k), 'timeLag:' num2str(l), ...
                          'coeff:' num2str(coeff)];
                           hillTerm(counter,:) = obj.weightMatrix(coeff ...
                               ,1) .* interval.^h/(1+interval.^h);
                       end
                       
                   end
               end
           end
           
            hillTerm = sum(hillTerm,1);
            figure;
            plot(interval, hillTerm, 'LineWidth', 1.5);
       end
   end
   
end