function convertAdjListGraphToPajekGraph(AdjGraph, filename)
% Converts an undirected but weighted adjacency list graph into an
% undirected Pajek graph, with node numbering starting from 1.
% Inputs:
%           AdjGraph -- n x 3 matrix (from_node, to_node, weight), node
%                       numbering can start either with 0 or 1.
%                       The format is assumed to be "two way", i.e. for each
%                       link between nodes a and b, there should be 
%                       (a,b,w) and (b,a,w) present in the list.
%                       The graph should contain no self loops.
%           filename -- string with filename of output file,
%                       e.g. 'pajekgraph.net'. If file already exists
%                       contents are overwritten
%
% last revision: 21/6/2011 by Michael  


% NOTE: the 0/1 functionaliy has not been thoroughly tested
% check if node numbering starts with 0 or 1
zero_or_one = min(AdjGraph(:,1));
if zero_or_one == 0
    display('Node numbering started with 0?!')
    % If node numbering starts with 0 increment
    AdjGraph(:,1) = AdjGraph(:,1) +1;
    AdjGraph(:,2) = AdjGraph(:,2) +1;
else
    display('Node numbering started with 1?!')
end

% find number of nodes
nr_nodes = max(AdjGraph(:,1));


% open file and write name list
fid = fopen(filename,'w+','n','ascii');
fprintf(fid,'*Vertices %i \n',nr_nodes);
for i= 1:nr_nodes
    fprintf(fid,'%i \"%i\" \n', int32(i), int32(i));
end

% delete redundant nodes
Index = AdjGraph(:,1)<AdjGraph(:,2);
Graph = AdjGraph(Index,:);

% start Edge section
fprintf(fid,'*Edges \n');

% write out pajek list..
for i=1:length(Graph)
    fprintf(fid,'%i %i %e \n', Graph(i,1), Graph(i,2), Graph(i,3));
end

fclose('all');


end
