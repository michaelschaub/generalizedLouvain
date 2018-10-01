function [] = stability_postprocess(mat_file, A,params)
%stability_postprocess -- postprocessing of stability results to obtain
%smoothed stability curves and best possible partitions from the obtained
%results.
% Input:    mat_file - path to matlab file in which stability results are
%                      stored eg 'matlab.mat' (created by save).
%           A        - adjacency matrix of the original graph
%           params   - struct storing the different options of stability
%requested variable names: S N C VI Time
%the effect of the function is to change the values in S, N, C_new and
%store into eg matlab_PP.mat (which can be recovered as
%load('matlab_PP.mat') )

% load stability files
load(mat_file);

if length(mat_file)<4 || ~strcmp(mat_file(end-3:end),'.mat')
    mat_file = [mat_file '.mat'];
end


% try to get PARAMS from data, if not possible and no input set to a default
if nargin<3 && ~exist('PARAMS','var')
    params.directed = false;
    params.teleport_tau = 0.15;
    params.laplacian = 'louvain_LNL';
    params.precomputed = false;
    params.NbNodes = length(A);
    warning('Could not infer parameters; default parameters have been set');
    
elseif exist('PARAMS','var') && nargin < 3
    params = PARAMS;
end


% copy old results and get set of all old found partitions
C_new = C;
C = unique(C','rows');
C = C';

% make numbering start from 1
if min(C)==0
    C= C+1;
end

% this definition needs to be here for transparency in parfor..
T = Time;
StabilityFunction = params.laplacian;


% Get transition matrices etc. if necessary
if ~params.precomputed
    params.precompute_only = true;
    [~, params] = feval(StabilityFunction, A, 0, params);
    params.precompute_only = false;
end
if params.ComputeParallel
    parfor j = 1:length(T)
        
        [solution, PARAMS] = feval(StabilityFunction, A, T(j), params);
        stabilities = zeros(1,size(C,2));
        
        for i=1:size(C,2)
            H = transformPartitionVectorToHMatrix(C(:,i));
            stabilities(i) = evaluate_clustering(H,solution,T(j),PARAMS);
        end
        
        % get best stability from the found partitions
        index = find(stabilities==max(stabilities),1);
        Stemp = stabilities(index);
        Ctemp = C(:,index);
        Ntemp = max(Ctemp);
        
        
        % store it in respective fields..
        S(j) = Stemp;
        C_new(:,j) = Ctemp;
        N(j) = Ntemp;
    end
    % store under old name, but append _PP.mat
    save([mat_file(1:end-4) '_PP.mat'],'S','N','VI','C_new','Time');
    
else % non parallel..
    for j=1:length(T)
        
        fprintf('Post-processing Markov time: %f\t(%i of %i)\n', T(j), ...
            j, length(Time))
        [solution, PARAMS] = feval(StabilityFunction, A, T(j), params);
        stabilities = zeros(1,size(C,2));
        
        for i=1:size(C,2)
            H = transformPartitionVectorToHMatrix(C(:,i));
            stabilities(i) = evaluate_clustering(H,solution,T(j),PARAMS);
        end
        
        % get best stability from the found partitions
        index = find(stabilities==max(stabilities),1);
        Stemp = stabilities(index);
        Ctemp = C(:,index);
        Ntemp = max(Ctemp);
        
        
        % store it in respective fields..
        S(j) = Stemp;
        C_new(:,j) = Ctemp;
        N(j) = Ntemp;
        % store under old name, but append _PP.mat
        save([mat_file(1:end-4) '_PP.mat'],'S','N','VI','C_new','Time');
    end
end % end else non parallel

end % end file


function stab = evaluate_clustering(H,solution,T,PARAMS)
if PARAMS.fixed_time
    
    if strcmp(PARAMS.louvain_type,'generalised')
        stab = trace(H' *solution*H);
        for kk = 1:2:size(PARAMS.null_vec,2)
            stab =  stab ...
                - (PARAMS.null_vec(:,kk)' *H)*(H'*PARAMS.null_vec(:,kk+1));
        end
        
    elseif strcmp(PARAMS.louvain_type,'normalised')
        stab = trace(H' *solution*H) ...
            - (PARAMS.pi' *H)*(H'*PARAMS.pi);
    end
    
elseif ~PARAMS.fixed_time
    
    if strcmp(PARAMS.louvain_type,'generalised')
        stab = trace(H' *T*solution*H);
        for kk = 1:2:size(PARAMS.null_vec,2)
            stab =  stab ...
                - (PARAMS.null_vec(:,kk)' *H)*(H'*PARAMS.null_vec(:,kk+1));
        end
        
    elseif strcmp(PARAMS.louvain_type,'normalised')
        stab = trace(H' *T*solution*H) ...
            - (PARAMS.pi' *H)*(H'*PARAMS.pi);
    end
    
    stab = (1-T) + stab;
    
else
    error('Something went wrong')
end

end
