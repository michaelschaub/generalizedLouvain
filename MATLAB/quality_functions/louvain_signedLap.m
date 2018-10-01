function [FlowMatrix, PARAMS] = louvain_signedLap(Graph, time, PARAMS)

PARAMS.louvain_type = 'generalised';
PARAMS.fixed_time = true;
% if PARAMS.directed == true
%     error('Not defined yet'); 
% end

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % store results for future use
    PARAMS.precomputed = true;   
    
    Dabs = diag(sum(abs(Graph),2));
    % signed Laplacian
    L = Dabs - Graph;
    
    if PARAMS.eig_decomp== false
        PARAMS.P = -L;
    elseif PARAMS.eig_decomp== true
        [V, lambda] = eig(full(-L')); % M' = V*lambda/V
        PARAMS.P = V;
        PARAMS.lambda = diag(lambda);
    else
        error('Something went wrong here')
    end
    
    
    if PARAMS.precompute_only
        FlowMatrix = 0;
        return
    end
    
end

N = length(Graph);
% "transition matrix" has been computed before
if PARAMS.eig_decomp == true
    % ATTENTION: exp(-L*t)^T = V*exp(lambda*t)*V^-1
    FlowMatrix = (PARAMS.P*diag(exp((PARAMS.lambda)*time))...
        /PARAMS.P);
    temp = real(FlowMatrix)*ones(N,1)/sqrt(N);
    PARAMS.null_vec = [temp, temp];
    
    FlowMatrix = real(FlowMatrix*FlowMatrix');
else
    % here P == L
    FlowMatrix = expm(time* PARAMS.P);
    temp = FlowMatrix'*ones(N,1)/sqrt(N);
    PARAMS.null_vec = [temp, temp];
    FlowMatrix = FlowMatrix'*FlowMatrix;
end
% symmetrization not needed but for numerical issues
FlowMatrix = (FlowMatrix +FlowMatrix')/2;
end
