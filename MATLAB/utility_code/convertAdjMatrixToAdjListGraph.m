function convertAdjMatrixToAdjListGraph(A,filename)
% Convert (sparse) adjacency matrix into adjacency list file that can be 
% used by the stability code. Node numbering is started with 0.
% Inputs:
% 
%           A:          Adjacency matrix of undirected graph 
%
%           filename:   string with filename of output file,
%                       e.g. 'stabilitygraph.txt'. If file already exists
%                       contents are overwritten
%
% last revision: 25/3/2010 by Michael  

fid = fopen(filename,'w+','n','ascii');

[i j link] = find(A);

for z=1:length(link)
    
    fprintf(fid,'%i %i %e \n', j(z)-1, i(z)-1, link(z));
end

fclose('all');

end
