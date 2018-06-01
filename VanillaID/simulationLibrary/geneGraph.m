classdef geneGraph
    properties (Access = public)
        simulationGraph
        numberOfGenes
        digraph
        highlight
        standardGroundTruth % Assumes 27 dictionary functions!!!
        adjacencyMatrix
    end
    methods
        function obj = geneGraph(numberOfGenes)
            obj.numberOfGenes = numberOfGenes;
            obj.adjacencyMatrix = zeros(numberOfGenes, numberOfGenes);
            obj.simulationGraph = simGraph(numberOfGenes);
            obj.digraph = digraph();
            for i=1:numberOfGenes
                obj.digraph = addnode(obj.digraph,{['Gene ' num2str(i)]});
            end
            obj.highlight = 0;
            obj.standardGroundTruth = zeros(27, numberOfGenes);
        end
        function obj = degradation(obj, whichGene, degradationRate)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) x, whichGene, whichGene, degradationRate);
            obj.standardGroundTruth(whichGene,whichGene) = degradationRate;
        end
        function obj = bias(obj, onWhichGene, biasGene, biasRate)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) 1, onWhichGene, biasGene, biasRate);
        end
        function obj = activation(obj, activated, activator, activationRate, hillOrder)
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) (x.^hillOrder)./(1+x.^hillOrder), activated, activator, activationRate);
            obj.digraph = addedge(obj.digraph,{['Gene ' num2str(activator)]}, {['Gene ' num2str(activated)]}, activationRate);
            obj.highlight = [obj.highlight 0];
            obj.standardGroundTruth(5*obj.numberOfGenes + ...,
                (hillOrder-1)*obj.numberOfGenes + activator,activated) = activationRate;
            obj.adjacencyMatrix(activator, activated) = activationRate;
        end
        function obj = repression(obj, repressed, repressor, repressionRate, hillOrder)
            % This is a simple repression without the option to add a K_M
            % constant
            obj.simulationGraph = obj.simulationGraph.addBasisFunction(@(x) 1./(1+x.^hillOrder), repressed, repressor, repressionRate);
            obj.digraph = addedge(obj.digraph,{['Gene ' num2str(repressor)]}, {['Gene ' num2str(repressed)]}, repressionRate);
            obj.highlight = [obj.highlight 1];
            obj.standardGroundTruth(obj.numberOfGenes + ...,
                (hillOrder-1)*obj.numberOfGenes + repressor,repressed) = repressionRate;
            obj.adjacencyMatrix(repressor, repressed) = -repressionRate;
            
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
        function printGraph(obj, fileName)
            fileID = fopen([fileName, '.tex'],'w');
            fprintf(fileID,'%s \n','\documentclass{article}');
            fprintf(fileID,'%s \n','\usepackage{tikz}');
            fprintf(fileID,'%s \n', '\usetikzlibrary{arrows,automata}');
            
            fprintf(fileID,'%s \n', '\begin{document}');
            fprintf(fileID,'%s \n', '\begin{tikzpicture}', ...,
                '[node distance=2.8cm, auto, line width = 0.5mm]');
            
            % Generate nodes
            fprintf(fileID,'%s (%s) {%s};\n', ...,
                '\node[shape=circle,draw=black]', ...,
                num2str(1), num2str(1));
            positionStrings = {'right of', 'below of', 'left of', 'below of'};
            
            for i=2:obj.numberOfGenes
                % If below, then previous line leader
                % If right of -> previous
                referenceString = num2str(i-1);
                positionIndex = mod(i+1,4);
                if (positionIndex == 0) 
                    positionIndex = 4;
                end
                fprintf(fileID,'%s (%s) [%s =%s] {%s};\n', ...,
                    '\node[shape=circle,draw=black]', ...,
                    num2str(i), positionStrings{positionIndex}, referenceString, ...,
                    num2str(i));
             
            end
            
            % Generate edges
            
            % Bend toggles
            bendToggle = ones(1,obj.numberOfGenes);
            bendStrings = {'right', 'left'};
            
            for i=1:obj.numberOfGenes
                for j=1:obj.numberOfGenes
                    if (obj.adjacencyMatrix(i,j) ~= 0)
                        currentBend = bendStrings{bendToggle(j)};
                        
                        % Consider activation first
                        if (obj.adjacencyMatrix(i,j) > 0)
                            
                            fprintf(fileID, ...,
                                '%s [->] (%s) edge[bend %s] node {$%s$} (%s);\n', ...,
                                '\path ', num2str(i), currentBend, ...,
                                num2str(obj.adjacencyMatrix(i,j)), num2str(j));
                        else
                            fprintf(fileID, ...,
                                '%s [-|] (%s) edge[bend %s] node {$%s$} (%s);\n', ...,
                                '\path ', num2str(i), currentBend,...,
                                num2str(obj.adjacencyMatrix(i,j)), num2str(j));
                        end
                        bendToggle(j) = bendToggle(j) + 1;
                        
                    end
                end
            end
            
            fprintf(fileID,'%s \n ','\end{tikzpicture}');
            
            fprintf(fileID,'%s \n','\end{document}');
            fclose(fileID);
            
            command = ['pdflatex ', fileName, '.tex'];
            [status,cmdout] = system(command);
            
            disp(cmdout);
        end
    end
end


