function showGraph(transitionMatrix)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here

   % % Find the number of nodes
   %  numNodes = size(transitionMatrix, 1);
   % 
   %  % Create dynamic node names using letters
   %  nodeNames = arrayfun(@(x) char('A' + (x - 1)), 1:numNodes, 'UniformOutput', false);
   % 
   %  % Find the source and target nodes for non-zero transitions
   %  [source, target] = find(transitionMatrix);
   % 
   %  % Define the edge weights (transition values)
   %  weights = nonzeros(transitionMatrix);
   % 
   %  % Create a directed graph from the source, target, and weights
   %  G = digraph(source, target, weights);
   % 
   %  % Assign dynamic node names
   %  G.Nodes.Name = nodeNames';
   % 
   %  % Plot the graph
   %  figure;
   %  h = plot(G, 'Layout', 'layered', 'NodeLabel', G.Nodes.Name);
   % 
   %  % Label the edges with the transition values
   %  labeledge(h, source, target, weights);
   % 
   %  % Add title and labels
   %  title('Graph Representation of Transitions');
   %  xlabel('Vertices');
   %  ylabel('Vertices');

  

% Dynamically create vertex names
numVertices = size(transitionMatrix, 1);
vertexNames = arrayfun(@(x) sprintf('V%d', x), 1:numVertices, 'UniformOutput', false);

% Create directed graph object
G = digraph(transitionMatrix, vertexNames);

% Plot the graph
figure;
h = plot(G, 'Layout', 'layered', 'EdgeLabel', G.Edges.Weight, ...
         'NodeFontSize', 12, 'EdgeFontSize', 10, 'ArrowSize', 15);

% Customize the graph appearance
title('State-Space Model / Automaton');
h.MarkerSize = 8; % Size of vertex markers
h.LineWidth = 1.5; % Thickness of edges
h.EdgeColor = '#4DBEEE'; % Edge color
h.NodeColor = 'blue'; % Node color
h.ArrowSize = 10; % Arrowhead size

grid on; % Add grid

end