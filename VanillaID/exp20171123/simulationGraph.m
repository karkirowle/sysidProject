classdef simulationGraph
   properties (Access = public)
       numberOfGenes
       hillOrder
       weightMatrix
   end
   methods
       function obj = simulationGraph(numberOfGenes,hillOrder)
           obj.numberOfGenes = numberOfGenes;
           obj.hillOrder = hillOrder;
           obj.weightMatrix = zeros(numberOfGenes + ...,
               hillOrder*2*numberOfGenes,numberOfGenes);
       end
       function obj = degradation(obj, degraded, amount)
           obj.weightMatrix(degraded,degraded) = -amount;
       end
       function obj = repression(obj, repressedBy, repressed, amount, order)
           % Total linear states + total previous order states + repressor
           % gene numbers
           row = obj.numberOfGenes ...
               + (order-1)*obj.numberOfGenes*2 ...
               + repressedBy;
           
           obj.weightMatrix(row,repressed) = amount;
       end
       function obj = activation(obj, activatedBy, activated, amount, order)
           % Total linear states + total previous order states + repressor
           % states + total activation numbers
           row = obj.numberOfGenes ...
               + (order-1)*obj.numberOfGenes*2 ...
               + obj.numberOfGenes ...
               + activatedBy;
           obj.weightMatrix(row,activated) = -amount;
       end
   end
   
end