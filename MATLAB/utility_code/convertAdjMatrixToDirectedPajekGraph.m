function convertAdjMatrixToDirectedPajekGraph(A,filename)
% Converts adjacency matrix of a network into pajek graph of an directed
% network. Assumes network to be fully connected.
% Inputs:
%           A:          Adjacency matrix of graph
%       
%           filename:   string with filename of output file,
%                       e.g. 'pajekgraph.net'. If file already exists
%                       contents are overwritten
%
% last revision: 30/8/2012 by Michael  

nr_nodes = length(A);
% open file and write name list
fid = fopen(filename,'w+','n','ascii');
fprintf(fid,'*Vertices %i \n',nr_nodes);
for i= 1:nr_nodes
    fprintf(fid,'%i \"%i\" \n', int32(i), int32(i));
end


% print arcs
fprintf(fid,'*Arcs\n');


[i j link] = find(A);
% writeout = [j i link];
% dlmwrite(filename,writeout,'-append','delimiter', ' ');

psize = length(link);
for z= 1:psize
            fprintf(fid,'%i %i %e \n', i(z), j(z), link(z));          

end
fclose('all');


end
