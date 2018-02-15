clc;
clear;
close all;

G = digraph();
e = G.Edges;
G = addnode(G,{'Gene 1'});
G = addnode(G,{'Gene 2'});
G = addnode(G,{'Gene 3'});
G = addedge(G,{'Gene 1'}, {'Gene 2'}, 10);
G = addedge(G,{'Gene 2'}, {'Gene 3'}, 10);
h = plot(G);
labeledge(h,1,{'10'})
highlight(h,1,2, 'EdgeColor', 'r')