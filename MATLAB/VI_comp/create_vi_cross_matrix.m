function [VImat] = create_vi_cross_matrix(dir_name)
% create vi matrix over time for a given directory containing all the partitions

cur_dir = pwd();

% augment dir_name if necessary
if ~strcmp(dir_name(1:11),'Partitions_')
    dir_name = ['Partitions_', dir_name];
end

cd(dir_name)
D = dir();
T = zeros(1,length(D));

for i=1:length(D)
    name = D(i).name;
    if length(name)>12 && strcmp(name(1:10),'Partition_')
        name_parts = strsplit(name,'_');
        timestr = name_parts{end};
        time = str2double(timestr(1:end-4));
        T(i) = time;
    end  
end

T= sort(T);
% there should be no zero in our times..
T(T==0) = [];
% ['Partition_' prefix '_' num2str(T(t),'%10.6f') '.mat']
VImat = zeros(length(T));
prefix = [strjoin({name_parts{1:end-1}},'_'), '_'];

for i = 1:length(T)
    name = [prefix, num2str(T(i),'%10.6f') '.mat'];
    comm1 = load(name);
    comm1 = comm1.communities;

    
    for j=i:length(T)
        name2 = [prefix, num2str(T(j),'%10.6f') '.mat'];
        comm2 = load(name2);
        comm2 = comm2.communities;
        vi_temp = varinfo_cross(comm1',comm2',true);
        VImat(i,j)= vi_temp;
    end
    
end
VImat = VImat+VImat'-diag(diag(VImat));
save('VIVsT','VImat');
cd(cur_dir);

end
