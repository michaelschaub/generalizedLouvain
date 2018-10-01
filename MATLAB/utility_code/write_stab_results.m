function write_stab_results(stab_file, outname, alt_ids)
% convenience function to quickly writeout the partitioning results of 
% stability in csv
if nargin <3
    alt_ids =[];
end

if nargin <2
    [~,name,~] = fileparts(stab_file);
    outname = [name '.csv'];
end

if length(outname)<4 || ~strcmp(outname(end-3:end),'.csv')
    outname = [outname '.csv'];
end

load(stab_file)

% check if we have a postprocessed file...
if(exist('C_new','var'))
    C = C_new;
end


n_times = length(Time);
prefixes = cell(1,n_times);
for i=1:n_times
    prefixes{i} = ['Time' num2str(Time(i))];
end

writeout_partition(C,outname,prefixes,alt_ids);

end