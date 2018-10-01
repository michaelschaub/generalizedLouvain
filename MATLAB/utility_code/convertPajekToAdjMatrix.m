function A= convertPajekToAdjMatrix(filename)
% Convert undirected pajek graph into sparse adjacency matrix. Assumes node
% numbering in pajek file starts with 0 or 1.
% Inputs:
%       
%           filename:   string with filename of output file,
%                       e.g. 'pajekgraph.net'. If file already exists
%                       contents are overwritten
% Outputs:
%           A:          Adjacency matrix of undirected graph in sparse data
%                       format 

% open file
fid = fopen(filename ,'r','n','ascii');

% skip "intro" part in pajek file
line=fgets(fid);
while(~feof(fid) && ~strcmp(line(1:6),'*Edges') ...
        && ~strcmp(line(1:6),'*edges') && ~strcmp(line(1:6),'*Arcs')...
        && ~strcmp(line(1:6),'*arcs'))  
    line=fgets(fid);
end

% scan adjacency list into G
G = fscanf(fid,'%e %e %e', [3 inf]);
fclose(fid);

% check if node numbering starts with 0 and in case make it start from 1
if min([G(1,:) G(2,:)]==0)
    G(1,:) = G(1,:)+1;
    G(2,:) = G(2,:)+1;
end

% find number of Nodes
N= max([G(1,:) G(2,:)]);

% allocate and fill in adjacency matrix
A= sparse(G(1,:),G(2,:),G(3,:),N,N);

% undirected pajek graphs are "one way"/nonsymmetric, hence symmetrize and remove double counting of self loops
A = A+A'-diag(diag(A));

end
