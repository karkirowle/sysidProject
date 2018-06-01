% Incoherent feed-forward loop

sim = geneGraph(3);

% Constructing basic interactions
sim = sim.activation(1, 2, 0.4, 4);
sim = sim.activation(2, 3, 0.4, 4);
sim = sim.repression(1, 3, 0.4, 4);

% Addding degradation terms
sim = sim.degradation(1, 0.1);
sim = sim.degradation(2, 0.1);
sim = sim.degradation(3, 0.1);

% Showing gene graph
sim.showGraph();

% Print graph
sim.printGraph('exampletikz');