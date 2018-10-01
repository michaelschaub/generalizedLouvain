function [vi,vi_mat] = varinfo_cross(partition_vectors1,partition_vectors2,ComputeParallel)
%compute the variation of information between two sets of partitions
%partition vectors correspond to C'
if nargin <3
    ComputeParallel = false;
end

number_of_partitions1 = size(partition_vectors1,1);    
number_of_partitions2 = size(partition_vectors2,1);
n = size(partition_vectors1,2);
n2 = size(partition_vectors2,2);

if n~=n2
    error('partitions should have same number of nodes')
end

vi_mat = zeros(number_of_partitions1,number_of_partitions2);
vi=0;

% check if partition vectors start with zero or one
if min(partition_vectors1(:))==0
    partition_vectors1 = partition_vectors1 + 1;
end
% check if partition vectors start with zero or one
if min(partition_vectors2(:))==0
    partition_vectors2 = partition_vectors2 + 1;
end 

% If all the partitions are identical, vi=0 and there is no need to do the
% rest of the calculations which are computationally expensive.
if  all(all([partition_vectors1; partition_vectors2]==...
        repmat(partition_vectors1(1,:),number_of_partitions1+number_of_partitions2,1)))
    return;
end

%Select only the partitions which are different 
[partition_vectors1,b,c1] = unique(partition_vectors1,'rows');
number_of_partitions1=length(b);
[partition_vectors2,b,c2] = unique(partition_vectors2,'rows');
number_of_partitions2=length(b);

vi_mat = zeros(number_of_partitions1,number_of_partitions2);

nodes = 1:n;

if nargin==3 && ComputeParallel
    parfor i = 1:number_of_partitions1
        partition_1 = partition_vectors1(i,:);
        partition_1 = double(partition_1);
        A_1 = sparse(partition_1,nodes,1);
        n_1_all = sum(A_1,2);
        vi_mat_row=vi_mat(i,:);

        for j = 1:number_of_partitions2
            partition_2 = partition_vectors2(j,:);
            partition_2 = double(partition_2);
            A_2 = sparse(nodes,partition_2,1);
            n_2_all = sum(A_2,1)';
            n_12_all = A_1*A_2;

            [rows,cols,n_12] = find(n_12_all);

            n_1 = n_1_all(rows);
            n_2 = n_2_all(cols);

            % fix this?!
            n_1 = reshape(n_1,size(n_12));
            n_2 = reshape(n_2,size(n_12));

            vi = sum(n_12.*log(n_12.^2./(n_1.*n_2)));
            vi = -1/(n*log(n))*vi;
            
            vi_mat_row(j)=vi;
            

        end
        vi_mat(i,:)=vi_mat_row;
    end
else
    for i = 1:number_of_partitions1
        partition_1 = partition_vectors1(i,:);
        partition_1 = double(partition_1);
        A_1 = sparse(partition_1,nodes,1);
        n_1_all = sum(A_1,2);

        for j = 1:number_of_partitions2
            partition_2 = partition_vectors2(j,:);
            partition_2 = double(partition_2);
            A_2 = sparse(nodes,partition_2,1);
            n_2_all = sum(A_2,1)';
            n_12_all = A_1*A_2;

            [rows,cols,n_12] = find(n_12_all);

            n_1 = n_1_all(rows);
            n_2 = n_2_all(cols);
            
            % fix this?!
            n_1 = reshape(n_1,size(n_12));
            n_2 = reshape(n_2,size(n_12));
            
            vi = sum(n_12.*log(n_12.^2./(n_1.*n_2)));
            vi = -1/(n*log(n))*vi;
            vi_mat(i,j)=vi;

        end
    end
end

vi_mat_new = zeros(number_of_partitions1,length(c2));

for i=1:number_of_partitions1
    vi_mat_new(i,:) = vi_mat(i,c2);
end
vi_mat=vi_mat_new(c1,:);

% vi_mat = vi_mat_new+vi_mat_new';
% vi_mat = vi_mat+vi_mat';

vi = mean(mean(vi_mat));

end
