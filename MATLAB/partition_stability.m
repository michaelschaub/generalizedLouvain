function [S, N, VI, C, PARAMS] = partition_stability(G, T, varargin)
%STABILITY    Graph partitioning optimizing stability with the Louvain
%             algorithm
%
%   [S, N, VI, C] = STABILITY(G, T) finds the optimal partitions of the 
%   graph G by optimizing the stability at each Markov time in vector T. G 
%   can either be the list of the edges in the graph (in the form [node i, 
%   node j, weight of link i-j; node k, node l, weight of link k-l;...] if 
%   the graph is weighted, or [node i, node j; node k, node l;...] if the 
%   graph is unweighted) or the adjacendy matrix of the graph. S, N, VI and 
%   C contain respectively the stability, the number of communities, the 
%   variation of information, and the optimal partition for each 
%   Markov time contained in T. If T is not specified, the modularity
%   (equivalent to stability for T=1) is calculated. Ideally, Markov time
%   should be sampled exponentially (e.g.: T = 10.^[-2:0.01:2]).
%
%   [S, N, VI, C] = STABILITY(G, T,'PARAM',VALUE) accepts one or more
%   comma-separated parameter name/value pairs. For a list of parameters 
%   and values, see "Parameter Options."
%
%    
%   Parameter Options:
% 
%        Parameter      Value                                 Default
%        ---------      -----                                 -------
%        L              Number of optimisations of the          100
%                       Louvain algorithm to be done at 
%                       each Markov time.
%
%        M              The top M partitions among the L        100
%                       given at each Markov time by the
%                       L louvain optimisations will be
%                       used to compute the variation of
%                       information.    
%
%        laplacian      Allows to choose which type of     'normalised/ louvain_FNL'
%                       quality function should be used.
%                       options include louvain_FNL, louvain_LNL,
%                       louvain_FCL, and louvain_LCL
%
%        directed       activate stability for directed         none
%                       graphs. Note that transition matrices
%                       are defined for left multiplications 
%                       here, i.e. A_ij is the link from i to j.
%
%        teleport_tau   teleportation probability               0.15
%                       (only active if directed == true)
%
%
%        noVI           Disables the calculation of the         none
%                       robustness of the partitions.
%                       Disabling this can significantly 
%                       speed up the calculations.
%
%        out            Enables saving step by step the         ''
%                       partitions found in a .mat file 
%                       located in the current folder, 
%                       along with the number of 
%                       communities (N), the value of 
%                       Stability (S), and the variation 
%                       of information (VI) for the 
%                       partitions obtained at each 
%                       Markov Time.
%
%       nocheck         Disables the checks for the             none
%                       encoding of the graph. This can 
%                       save computational time but can 
%                       also lead to serious errors if 
%                       the graph has not been properly
%                       encoded.
%
%       prec            Precision: defines a threshold for      1e-18
%                       the range of weights allowed in  
%                       the laplacian exponential matrix 
%                       of the full stability.
%
%       plot            Plots the plots the results of          none
%                       the stability, number of 
%                       communities and variation of 
%                       information as a function of the
%                       Markov time.           
%
%       v               Verbose mode.                           none
%
%       p               Parallel mode.                          none
%
%       t                                                       none
%                       Enables saving step by step the         
%                       partitions found at each Markov 
%                       time in a file located in a 
%                       folder named 'Partitions_..' 
%                       The option 'out' must be on.
%
%


% Unparsed default parameters
Time = 1;                                       % Markov times at which the graph should be partitioned
flag_matlabpool = false;                        % for opening/closing workpool for parallel computation

if nargin<1    
    error(['Please provide the graph to be partitioned. Type "help'...
        'stability" for more information.']);
end

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$          Arguments parsing               $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
[Graph, StabilityFunction, OutputFile, prefix, Sanity, plotStability, verbose, TextOutput, PARAMS]... 
    = parseinput(length(varargin),G,varargin);

%----------------------
% Check Argument 2: T
if nargin > 1
    if (isvector(T) && isnumeric(T))
        Time=T;
    else
        error(['The second argument should be a numerical vector.'...
            'Type "help stability" for more information.']);
    end
end

% Initialisation
S = zeros(1, length(Time));
N = zeros(1, length(Time));
VI = zeros(1, length(Time));
C = zeros(PARAMS.NbNodes, length(Time));

% Parallelization Initialize if matlabpool is not yet running.

if PARAMS.ComputeParallel && isempty(gcp)
    flag_matlabpool = true;
    parpool
end

if TextOutput
    mkdir(['Partitions_' prefix]);
end
if OutputFile
    save(['Stability_' prefix '.mat'],'Time','S','N','VI','C','PARAMS');
end
if plotStability
    figure_handle = figure;
end
if verbose
    step_prec=0;
end

%-------------------------
% End Initialization
%-------------------------



% Print Information
if verbose
    disp(' ');
    disp('   Stability will be computed with the following parameters:');
    disp(' ');
    disp(['      - ' StabilityFunction]);
    if PARAMS.directed 
        disp('      - DIRECTED graph'); 
        disp(['      - teleportation parameter set to ' num2str(PARAMS.teleport_tau) '(may be unused)']); 
    end
    if PARAMS.ComputeVI; disp('      - Computation of the variation of information (VI): Yes'); else disp('      - Computation of the variation of information: No'); end
    if OutputFile; disp(['      - Save the results in a file: Yes ' prefix]); else disp('      - Save the results in a file: No'); end
    if OutputFile; disp(['      - Saved file prefix: ' prefix]); end
    if Sanity; disp('      - Check the input graph: Yes'); else disp('      - Check the input graph: No'); end
    if plotStability; disp('      - Plot the results: Yes'); else disp('      - Plot the results: No'); end
    if verbose; disp('      - Verbose mode: Yes');    else disp('      - Verbose mode: No');    end
    if PARAMS.ComputeParallel; disp('      - Parallel computation: Yes'); else disp('      - Parallel computation: No'); end
    disp(['      - Number of Louvain iterations: ' int2str(PARAMS.NbLouvain)]);
    disp(['      - Number of Louvain iterations used for the computation of VI: ' int2str(PARAMS.M)]);
    disp(['      - Weight Precision used: ' num2str(PARAMS.Precision)]); 
    disp(['      - Louvain Precision used: ' num2str(PARAMS.PrecisionLouvain)]); 
    disp(' ');
    tstart=tic;
end





%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$      Computation of  stability           $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

% Get transition matrices etc.
if ~PARAMS.precomputed
    PARAMS.precompute_only = true;
    [~, PARAMS] = feval(StabilityFunction, Graph, 0, PARAMS);
    PARAMS.precompute_only = false;
end

% Loop over all Markov times
if (PARAMS.ComputeParallel && ~PARAMS.eig_decomp)
    if TextOutput
        cd(['Partitions_' prefix]);
    end
    parfor t=1:length(Time)        
        if verbose
            disp(['   Partitioning for Markov time = '...
                num2str(Time(t),'%10.6f') '...']);
        end

        [FlowMatrix, tempP] = feval(StabilityFunction, Graph, Time(t), PARAMS);

        %prune out weights that are too small as defined by precision
        FlowMatrix=max(max(FlowMatrix))*PARAMS.Precision*...
            round(FlowMatrix/(max(max(FlowMatrix))*PARAMS.Precision));
        
        %change into list format for optimisation with Louvain
        [row,col,val] = find(tril(FlowMatrix));
        graph=[col-1,row-1,val];
        
        
        %Optimize with Louvain NbLouvain times
        if PARAMS.fixed_time
            mtime =1;
        else 
            mtime = Time(t);
        end
        if isfield(tempP,'null_vec')
            null_vec = tempP.null_vec;
                    else 
            null_vec = [];
        end
        [stability, nb_comm, communities] = ...
            louvain(graph, mtime, PARAMS.NbLouvain, ...
            PARAMS.PrecisionLouvain,PARAMS.louvain_type,randi(intmax),...
            null_vec);
        
        
        %Comment: maybe one should pick one of the best solutions at random,
        %if two solutions have the same value;
        index = find(stability==max(stability),1);
        
        S(t) = stability(index);
        C(:,t) = communities(:,index);
        N(t) = nb_comm(index);
        
        if PARAMS.ComputeVI
            VI(t) = computeRobustness(communities, stability, PARAMS.M,PARAMS.ComputeParallel);
        else
            VI(t)=0;
        end
        
        
        if TextOutput
            communities = int32(communities);
            parsave(['Partition_' prefix '_' num2str(Time(t),'%10.6f') '.mat'],communities);
        end
        
    end
    
    if plotStability
        stability_plot(Time,Time(end),S,N,VI,PARAMS.ComputeVI,figure_handle);
    end
    
    if TextOutput
        cd ..
    end
        
    if OutputFile
        save(['Stability_' prefix '.mat'],'PARAMS','Time','S','N','VI','C','-append');
    end
    
    
else 
    
    for t=1:length(Time)
        
        if verbose
            disp(['   Partitioning for Markov time = '...
                num2str(Time(t),'%10.6f') '...']);
        end
        
        
        [FlowMatrix, PARAMS] = feval(StabilityFunction, Graph, Time(t), PARAMS);
        
        
        % prune out weights that are too small as defined by precision
        FlowMatrix=max(max(FlowMatrix))*PARAMS.Precision*...
            round(FlowMatrix/(max(max(FlowMatrix))*PARAMS.Precision));
        
        % change into list format for optimisation with Louvain
        [row,col,val] = find(tril(FlowMatrix));
        clear FlowMatrix
        graph=[col-1,row-1,val];
        
        if PARAMS.fixed_time
            mtime =1;
        else
            mtime = Time(t);
        end
        if isfield(PARAMS,'null_vec')
            null_vec = PARAMS.null_vec;
        else 
            null_vec = [];
        end

        
        % Optimize with Louvain NbLouvain times
        if PARAMS.ComputeParallel && PARAMS.eig_decomp
            nr_threads = parpool('size');
            stability = cell(1,nr_threads);
            nb_comm_temp = cell(1,nr_threads);
            communities = cell(1,nr_threads);
            shares = split_even(PARAMS.NbLouvain,nr_threads);
            % computation in parallel with cell arrays
            parfor l=1:nr_threads;
                [stability{l}, nb_comm_temp{l}, communities{l}] = ...
                    louvain(graph, mtime, shares(l), PARAMS.PrecisionLouvain ,...
                    PARAMS.louvain_type ,randi(intmax),null_vec);
            end
            % assignements
            communities = cat(2,communities{:});
            stability = cat(2,stability{:});
            nb_comm = cat(2,nb_comm_temp{:});
            
        % non parallel version
        else
            % Optimize with Louvain NbLouvain times
            [stability, nb_comm, communities] = ...
                louvain(graph, mtime, PARAMS.NbLouvain, ...
                PARAMS.PrecisionLouvain,PARAMS.louvain_type,randi(intmax),...
                null_vec);
        end 
        clear graph;
        
        % Comment: maybe one should pick one of the best solutions at random,
        % if two solutions have the same value;
        index = find(stability==max(stability),1);
        
        S(t) = stability(index);
        C(:,t) = communities(:,index);
        N(t) = nb_comm(index);
        if PARAMS.ComputeVI
            VI(t) = computeRobustness(communities, stability, PARAMS.M,PARAMS.ComputeParallel);
        else
            VI(t)=0;
        end
        
        
        if plotStability && t>1
            stability_plot(Time,t,S,N,VI,PARAMS.ComputeVI,figure_handle);
        end
        
        if TextOutput
            cd(['Partitions_' prefix]);
            comm = int32(communities);
            save(['Partition_' prefix '_' num2str(Time(t),'%10.6f') '.mat'],'comm');
            cd ..;
        end
        
        if OutputFile
            save(['Stability_' prefix '.mat'],'PARAMS','Time','S','N','VI','C','-append');
        end
        
        
        if verbose && 100*t/length(Time) >= step_prec+10
            disp(' ');
            disp(['   Completed: ' num2str(round(100*t/length(Time)),10) '%']);
            remaining_time=toc(tstart)*(1-t/length(Time))/(t/length(Time));
            nb_hours = floor(remaining_time/3600);
            nb_min = floor((remaining_time - nb_hours*3600)/60);
            nb_sec = round(remaining_time - nb_hours*3600 - nb_min*60);
            disp(['   Estimated time remaining: ' datestr([2011  1 1 nb_hours nb_min nb_sec], 'HH:MM:SS')]);
            disp(' ');
            step_prec = floor(100*t/length(Time));
        end
        
    end
end


if verbose
    c = clock;
    disp(' ');
    disp(['   Partitioning of the graph finished at ' datestr([2011 1 1 c(4) c(5) c(6)], 'HH:MM:SS')]);
    remaining_time=toc(tstart);
    nb_hours = floor(remaining_time/3600);
    nb_min = floor((remaining_time - nb_hours*3600)/60);
    nb_sec = round(remaining_time - nb_hours*3600 - nb_min*60);
    disp(['   Total time needed: ' datestr([2011 1 1 nb_hours nb_min nb_sec], 'HH:MM:SS')]);
end

if flag_matlabpool
    parpool close;
end



end


%------------------------------------------------------------------------------
% HELPER FUNCTION: PARSE ARGUMENT
%------------------------------------------------------------------------------
function [G, StabilityFunction, OutputFile, prefix, Sanity, plotStability, verbose, TextOutput, PARAMS]...
        = parseinput(options,G,varargin)
% Parse the options from the command line

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$          Default parameters              $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Global" options relevant for output and control flow

stability_type_specified=false;
threshold_nnodes_full = 1000;
threshold_nedges_full = 5000;

OutputFile = false;                             % No output file by default.

PARAMS.laplacian = 'louvain_FNL';               % Default Laplacian dynamics
Sanity = true;                                  % If true, performs the graph sanity checks
plotStability = false;                          % If true, plots the results of the stability, 
                                                % number of communities and variation of 
                                                % information vs Markov time.           
verbose = false;                                % Toggle verbose mode
prefix = '';                                    % Output prefix
TextOutput = false;                             % Toggles the text output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options stored in struct relevant for optimization etc.
PARAMS = struct;                                % create empty structure for storing parameters
PARAMS.precomputed = false;                     % Flag for precomputed transition matrix + stationary distribution
PARAMS.precompute_only = false;                 % set to true if you only want to compute the transition matrix + stat. dist
PARAMS.directed = false;                        % enables dealing with directed graphs
PARAMS.ComputeVI = true;                        % True if the variation of information should be computed
PARAMS.ComputeES = false;                       % True if edge statistics should be computed
PARAMS.ComputeParallel = false;                 % Toggles the computation in parallel
PARAMS.NbLouvain = 100;                         % Number of louvain optimisations at each Markov time
PARAMS.NbNodes = 0;                             % Total number of nodes;
PARAMS.Precision = 1e-18;                       % Threshold edges weigths
PARAMS.PrecisionLouvain = 1e-12;                % Threshold for minimal Louvain improvement
PARAMS.M = 100;                                 % Top M partitions among the L found by louvain 
                                                % are used to compute the variation of information
PARAMS.K = NaN;                                 % K stabilities value, only relevant for Ruelle random walk 
                                                % and k stabilities. K =-1 corresponds to the normalised Laplacian
PARAMS.teleport_tau = 0.15;                     % teleportation probability (only relevant for directed graphs)
PARAMS.eig_decomp = false;%true;                % use eigendecomposition or not

    

% --------------------------
% actual parsing begins here
attributes={'novi', 'l', 'm', 'out', 'nocheck', 'laplacian', 'prec','precLouvain', 'plot','v','t','p','k','directed','teleport','eig'};

if options > 0
    
    varargin = varargin{:}; % extract cell array input from varargin
    
    % test whether attribute-value pairs are specified, or fixed parameter order
    stringoptions = lower(varargin(cellfun('isclass',varargin,'char')));
    attributeindexesinoptionlist = ismember(stringoptions,attributes);
    newinputform = any(attributeindexesinoptionlist);
    if newinputform
        % parse values to functions parameters
        i = 1;
        while (i <= length(varargin))
            if strcmpi(varargin{i},'directed')
                PARAMS.directed = true;
                i = i+1;
            elseif strcmpi(varargin{i},'novi')
                PARAMS.ComputeVI = false;
                i = i+1;
            elseif strcmpi(varargin{i},'nocheck')
                Sanity = false;
                i = i+1;
            elseif strcmpi(varargin{i},'plot')
                plotStability = true;
                i = i+1;
            elseif strcmpi(varargin{i},'v')
                verbose = true;
                i = i+1;
            elseif strcmpi(varargin{i},'eig')
                PARAMS.eig_decomp = true;
                i = i+1;
            elseif strcmpi(varargin{i},'p')
                if exist('parpool','file')
                    PARAMS.ComputeParallel = true;
                else
                    PARAMS.ComputeParallel = false;
                    warning('The Parallel Computing Toolbox of Matlab does not appear to be installed. Defaulting to single node computation...');
                end
                i = i+1;
            elseif strcmpi(varargin{i},'t')
                TextOutput = true;
                i = i+1;
            else
                %Check to make sure that there is a pair to go with
                %this argument.
                if length(varargin) < i + 1
                    error('MATLAB:stability:AttributeList', ...
                        'Attribute %s requires a matching value', varargin{i});
                elseif strcmpi(varargin{i},'laplacian')
                    if ischar(varargin{i+1})
                        Laplacian = varargin{i+1};
                        stability_type_specified = true;
                    else
                        error('MATLAB:stability:laplacian',...
                            'Please provide a matching value for attribute laplacian. It must either be ''normalised'' or ''combinatorial''.');
                    end
                elseif strcmpi(varargin{i},'l')
                    if isnumeric(varargin{i+1})
                        PARAMS.NbLouvain = round(varargin{i+1});
                        PARAMS.M = round(varargin{i+1});
                    end
                elseif strcmpi(varargin{i},'prec')
                    if isnumeric(varargin{i+1})
                        PARAMS.Precision = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'precLouvain')
                    if isnumeric(varargin{i+1})
                        PARAMS.PrecisionLouvain = varargin{i+1};
                    end               
                elseif strcmpi(varargin{i},'m')
                    if isnumeric(varargin{i+1})
                        PARAMS.M = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'k')
                    if isnumeric(varargin{i+1})
                        PARAMS.K = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'teleport')
                    if isnumeric(varargin{i+1})
                        PARAMS.teleport_tau = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'out')
                    if ischar(varargin{i+1})
                        OutputFile = true;
                        prefix = varargin{i+1};
                    else
                        error('MATLAB:stability:out',...
                            'Please provide a matching value for attribute out. It must be a string.');
                    end
                else
                    error('MATLAB:stability:Attribute',...
                        'Invalid attribute tag: %s', varargin{i});
                end
                i = i+2;
            end
        end
    else 
        if ischar(varargin{1})
            error('MATLAB:stability:Attribute',...
                            'Invalid attribute tag: %s', varargin{1});
        else
            error('MATLAB:stability:Attribute',...
                            'Invalid attribute tag: %d', varargin{1});
        end
    end
end

TextOutput = TextOutput && OutputFile;

if ~stability_type_specified
        disp('  ---------------------------------------------------------------');
        disp('  You have not specified the quality function / dynamics ');
    if (size(G,1) == size(G,2) && nnz(G)<threshold_nedges_full && size(G,1)<threshold_nnodes_full) ...
            || (size(G,1) ~= size(G,2) && size(G,1)<threshold_nedges_full && max(max(G(:,1:2)))<threshold_nnodes_full)
        Laplacian = 'louvain_FNL';
        disp('  Stability with normalised Laplacian dynamics will be used');
        disp('  ---------------------------------------------------------------');
    else
        Laplacian = 'louvain_LNL';
        disp('  Linearised stability with normalised Laplacian dynamics will be used');
        disp('  ---------------------------------------------------------------');
    end
 
end

% input defined which stability function is used
if exist(Laplacian,'file')==2
    StabilityFunction = Laplacian;
else
    error('Please provide a valid matching value for attribute laplacian.');
end

PARAMS.laplacian = Laplacian;

%-----------------------
% Argument 1: G
    
% Transform graph input matrix representation
% first case -- matrix input
if size(G,1) == size(G,2) 
    Graph=G;
    % second case -- edgelist input
elseif size(G,1) ~= size(G,2)
    
    % check if node numbering starts with zero or one
    % if zero we need to add 1, to make it a proper matlab index
    min_add = 1 - min([G(:,1);G(:,2)]);
    num_nodes = max([G(:,1);G(:,2)])+min_add;
    
    if size(G,2)==3
        Graph=sparse(G(:,1)+min_add,G(:,2)+min_add,G(:,3),...
            num_nodes,num_nodes);
    elseif size(G,2)==2
        Graph=sparse(G(:,1)+min_add,G(:,2)+min_add,1,...
            num_nodes,num_nodes);
    else
        error(['Wrong size for G: G should be a graph '...
            'saved either as a list of edges (size(G)=[N,3] '...
            'if weighted, size(G)=[N,2] if unweighted) or as an ' ...
            'adjacency matrix (size(G)=[N,N])']);
    end
    
end

% get number of nodes and edges
PARAMS.NbNodes = size(Graph,2);

% Check if the graph is correctly encoded
if Sanity
    G=check(Graph, verbose, PARAMS);
end

end



%------------------------------------------------------------------------------
% HELPER FUNCTION TO DISTRIBUTE LOAD ACROSS MULTIPE THREADS
%------------------------------------------------------------------------------
function shares = split_even(N,nr_threads)
%Function to compute an even split of N runs between t threads 
shares = ones(1,nr_threads)*floor(N/nr_threads);
shares(1) = shares(1)+ rem(N,nr_threads);
end

%------------------------------------------------------------------------------
% HELPER FUNCTION TO SAVE INSIDE PARFOR
%------------------------------------------------------------------------------
function parsave(varargin)

savefile = varargin{1}; % first input argument
savevar = struct;
for i=2:nargin
    savevar.(inputname(i)) = varargin{i}; % other input arguments
end
save(savefile,'-struct','savevar')

end


%------------------------------------------------------------------------------
% HELPER FUNCTION TO CHECK WHETHER INPUT GRAPH IS ACCEPTABLE
%------------------------------------------------------------------------------
function Graph = check(Graph, verbose, PARAMS)
% Check that the graph is properly encoded.
    if verbose
        disp(' ');
        disp('   Graph sanity check...');
    end
        
    if size(Graph,2) ~= size(Graph,1)
            error('Graph has wrong size');
    end

    % Check that graph contains just numbers
    if any(any(~isnumeric(Graph)))
        error(['The graph provided contains elements which are not' ... 
            'numbers (isnumeric == false). Please check your graph' ... 
            'and try again.']);
    end
        
    % Check symmetry of the adjacency matrix if graph is not directed
    if PARAMS.directed == false
        if any(any(Graph~=Graph'))
            if nnz(triu(Graph,1))>0 && nnz(tril(Graph,-1))>0
                error('The graph provided is a directed graph.');
            else
                warning(['Adjacency matrix A of provided graph is'...
                    'triangular -- symmetrizing A = A + A^T']);
                Graph=Graph+Graph';
            end
        end
    end
    
    % Check for isolated nodes
    if ( any( sum(abs(Graph))' == 0 & sum(abs(Graph),2) == 0 ) )
        warning('There are isolated nodes in the graph!?');
    end
    
    % Check for disconnected components
    if exist('graphconncomp','file') == 2
        nbcomp=graphconncomp(sparse(Graph),'WEAK',true);
        if nbcomp>1
            warning(['There are ' num2str(nbcomp) '(weakly) connected '...
                'components in the graph. If your graph is directed '...
                'please be aware of the teleportation settings.']);
        end
    end
end


%------------------------------------------------------------------------------
% HELPER FUNCTION TO COMPUTE ROBUSTNESS VALUES
%------------------------------------------------------------------------------
function VI = computeRobustness(lnk, lnkS, M,ComputeParallel)

% Parameter
[~,i] = sort(lnkS);
lnk=lnk(:,i);
lnk=lnk(:,end-M+1:end);
VI = varinfo(lnk',ComputeParallel);
clear i;
end


%------------------------------------------------------------------------------
% HELPER FUNCTION FOR PLOTTING ON THE FLY
%------------------------------------------------------------------------------
function [] = stability_plot(Time,t,S,N,VI,ComputeVI,figure_handle)

set(0,'CurrentFigure',figure_handle);

if ComputeVI
    subplot(2,1,1), ax=plotyy(Time(1:t),N(1:t),Time(N>1),S(N>1));
else
    ax=plotyy(Time(1:t),N(1:t),Time(N>1),S(N>1));
end
xlabel('Markov time');
set(ax(1),'YScale','log');
set(ax(2),'YScale','log');
set(ax(1),'YTickMode','auto','YTickLabelMode','auto','YMinorGrid','on');
set(ax(2),'YTickMode','auto','YTickLabelMode','auto','YMinorGrid','on');
set(get(ax(1),'Ylabel'),'String','Number of communities');
set(get(ax(2),'Ylabel'),'String','Stability');
set(ax(1),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [1 10^ceil(log10(max(N)))], 'XScale','log','XMinorGrid','on');
set(ax(2),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [10^floor(log10(min(S(N>1)))), 1], 'XScale','log');
ylabel('Number of communities');
if ComputeVI 
    subplot(2,1,2), semilogx(Time(1:t),VI(1:t));
    set(gca, 'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YMinorGrid','on','XMinorGrid','on');
    if max(VI)>0
        set(gca,'YLim', [0 max(VI)*1.1]);
    end
    xlabel('Markov time');
    ylabel('Variation of information');
end
drawnow;

end


