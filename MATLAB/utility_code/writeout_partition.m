function writeout_partition(C,filename,prefixes,alternative_ids)
%function for creating a csv file assigning the egdes to communities
% C: clustering matrix over time
% filename: file for write-out
% prefixes: cell array containing the names of the categories, if not given
% time-stamps are used
% alternative ids: If original graph numbering was different provide a map
% to the alternative ids as row vector

[nr_nodes, nr_time] = size(C);


if nargin<3 || isempty(prefixes)
    t_indices = 1:nr_time;
    % put headers of columns
    fid = fopen(filename,'w+');
    fprintf(fid,'Id');
    fprintf(fid,' time%1d',t_indices);
    fprintf(fid,'\n');
    fclose(fid);
else
    HEADER = horzcat('Id', prefixes);
    dlmcell(filename,HEADER,' ');
end

DATA = cell(nr_nodes,nr_time+1);

if nargin <4 || isempty(alternative_ids)
    DATA(:,1) = num2cell((1:nr_nodes)');
else
    DATA(:,1) = num2cell(alternative_ids');
end

% type in column 3 + clusterings
DATA(:,2:end) = num2cell(C);

%write stuff in file
dlmcell(filename,DATA,' ','-a');



end
