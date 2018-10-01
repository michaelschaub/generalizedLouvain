function [vi,vi_mat] = varinfo(partition_vectors,ComputeParallel)
%VARINFO      Calculates the variation of information matrix and average
%             between all pairs of a set of partitions
%
%   [VI,VI_MAT] = VARINFO(P) calculates the variation of information between
%   each pair of partitions contained in P, where P is the M by N matrix
%   of partitions where N is the number of nodes in the original graph
%   and M is the number of partitions. The output VI is the average variation
%   of information between all pairs of partitions, and VI_MAT is the M by M 
%   matrix where entry (i,j) is the variation of information between 
%   the partitions contained in column i and j of the matrix P.
%   
%   [VI,VI_MAT] = VARINFO(P,F) allows the calculation of the variation of
%   information in parallel if the boolean F is true and provided that
%   matlab pool is running.
%
%   This code has been adapted from the code originally implemented for the 
%   following paper:
%
%       The performance of modularity maximization in practical contexts.
%       B. H. Good, Y.-A. de Montjoye and A. Clauset.
%       Physical Review E 81, 046106 (2010).
%
%   The original code can be found at: 
%   http://tuvalu.santafe.edu/~aaronc/modularity/

if nargin <2
    ComputeParallel = false;
end

number_of_partitions = size(partition_vectors,1);    
n = size(partition_vectors,2);

vi_mat = zeros(number_of_partitions);
vi=0;
 
% check if partition vectors start with zero or one
if min(partition_vectors(:))==0
    partition_vectors = partition_vectors + 1;
end    
% If all the partitions are identical, vi=0 and there is no need to do the
% rest of the calculations which are computationally expensive.
if  all(all(partition_vectors==repmat(partition_vectors(1,:),number_of_partitions,1)))
    return;
end

%Select only the partitions which are different 
[partition_vectors,b,c] = unique(partition_vectors,'rows');

number_of_partitions=length(b);

vi_mat = zeros(number_of_partitions);

nodes = 1:n;

if nargin==2 && ComputeParallel
    parfor i = 1:number_of_partitions
        partition_1 = partition_vectors(i,:);
        partition_1 = double(partition_1);
        A_1 = sparse(partition_1,nodes,1);
        n_1_all = sum(A_1,2);
        vi_mat_row=vi_mat(i,:);

        for j = 1:i-1
            partition_2 = partition_vectors(j,:);
            partition_2 = double(partition_2);
            A_2 = sparse(nodes,partition_2,1);
            n_2_all = sum(A_2,1)';
            % this is the (unnormalised) confusion matrix!
            n_12_all = A_1*A_2;
            
            temp = diag(n_1_all) \ (n_12_all.^2) / diag(n_2_all);

%             logterm = spfun(@log,temp);
            vi = -1/(n*log(n))*sum(nonzeros(n_12_all).*log(nonzeros(temp)));           
            
            vi_mat_row(j)=vi;
            
        end
        vi_mat(i,:)=vi_mat_row;
    end
else
    for i = 1:number_of_partitions
        partition_1 = partition_vectors(i,:);
        partition_1 = double(partition_1);
        A_1 = sparse(partition_1,nodes,1);
        n_1_all = sum(A_1,2);

        for j = 1:i-1
            partition_2 = partition_vectors(j,:);
            partition_2 = double(partition_2);
            A_2 = sparse(nodes,partition_2,1);
            n_2_all = sum(A_2,1)';
            n_12_all = A_1*A_2;

            temp = diag(n_1_all) \ (n_12_all.^2) / diag(n_2_all);

            vi = -1/(n*log(n))*sum(nonzeros(n_12_all).*log(nonzeros(temp)));
            
            vi_mat(i,j)=vi;

        end
    end
end

vi_mat_new = zeros(number_of_partitions,length(c));

for i=1:number_of_partitions
    vi_mat_new(i,:) = vi_mat(i,c);
end
vi_mat_new=vi_mat_new(c,:);

vi_mat = vi_mat_new+vi_mat_new';

vi = mean(squareform(vi_mat));

end
