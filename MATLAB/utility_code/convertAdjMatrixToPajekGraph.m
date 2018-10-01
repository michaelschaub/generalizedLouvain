function convertAdjMatrixToPajekGraph(A,filename)
% Converts adjacency matrix of undirected network into pajek graph of an
% undirected network. Assumes network to be fully connected.
% Inputs:
%           A:          Adjacency matrix of undirected graph
%       
%           filename:   string with filename of output file,
%                       e.g. 'pajekgraph.net'. If file already exists
%                       contents are overwritten
%
% last revision: 25/3/2010 by Michael  

nr_nodes = length(A);
% open file and write name list
fid = fopen(filename,'w+','n','ascii');
fprintf(fid,'*Vertices %i \n',nr_nodes);
for i= 1:nr_nodes
    fprintf(fid,'%i \"%i\" \n', int32(i), int32(i));
end


% print edges
fprintf(fid,'*Edges\n');
% This is important due to column major storage of Matlab
A = tril(A);
[i j link] = find(A);
psize = length(link);
percent = 10;
for z= 1:psize
            % looping over lower triangular hence swapping indices, but A=A'!
            fprintf(fid,'%i %i %e \n', j(z), i(z), link(z));
            
            if z >= psize*percent/100
                fprintf('%i percent done\n',percent)
                percent = percent+10;
            end
end

fclose('all');

end
