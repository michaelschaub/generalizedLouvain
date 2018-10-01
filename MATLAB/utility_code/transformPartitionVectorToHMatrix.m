function H = transformPartitionVectorToHMatrix(pvector)
%Transforms a given partition vector into a valid H matrix (Partition incidence matrix)
% assumes that the ordering is contigous and starts from either zero or
% one.

if min(pvector)==0
    pvector = pvector+1;
end
nr_nodes = length(pvector);

% create H matrix
H = sparse(1:nr_nodes,pvector,1);


end
