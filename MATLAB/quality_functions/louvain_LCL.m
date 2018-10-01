function [FlowMatrix, PARAMS] = louvain_LCL(Graph, time, PARAMS)
% Computes the linearised combinatorial stabilty matrix
%
PARAMS.louvain_type = 'generalised';

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        PARAMS.fixed_time = true;
        dout = sum(Graph,2);
        av_deg = mean(dout);
        dangling = (dout==0);
        dout(dangling) = 1;
        Dout = sparse(diag(dout));
        clear dout;
        M1 = (Dout./av_deg)*(1-PARAMS.teleport_tau)*(Dout\Graph); 
        % deterministic part of transition
        
        a = (Dout./av_deg)*...
            diag(PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes,1);
        b = ones(PARAMS.NbNodes,1)/PARAMS.NbNodes;

        
        M = (M1 + a*b');
        clear dangling
        L = ( Dout/av_deg ) - M;
        % BEWARE, this might be slow!
        [V, lambda] = eigs(L',5,'SR'); % L' = V*lambda/V
        lambda_max = min(diag(real(lambda)));
        v = V(:,diag(lambda) == lambda_max);
        v = abs(v);              % make sure eigenvector is positive
        % store results for future use
        PARAMS.precomputed = true;
        PARAMS.pi = v/sum(v);
        PARAMS.P = diag(PARAMS.pi)*(M1 -Dout/av_deg);
        
        PARAMS.u = -diag(PARAMS.pi)*a;
        PARAMS.v = b;
        
        if PARAMS.precompute_only
            FlowMatrix = 0;
            return
        end
        
        
        
    else % undirected case
        
        PARAMS.fixed_time = false;
        
        dout = sum(Graph,2);
        av_deg = mean(dout);

        % store results for future use
        PARAMS.precomputed = true;       
        % stat. dist is uniform
        PI= sparse(1:PARAMS.NbNodes,1:PARAMS.NbNodes,1/PARAMS.NbNodes);  
        PARAMS.pi = diag(PI);   % store results for future use      
        PARAMS.P = Graph./av_deg;    
        
        if PARAMS.precompute_only
            FlowMatrix = 0;
            return
        end
        
    end
    
end

% stationary distribution and "transition matrix" have been computed before
if PARAMS.directed == true
    FlowMatrix = diag(PARAMS.pi) + PARAMS.P*time;
    PARAMS.null_vec = [PARAMS.pi, PARAMS.pi, time*PARAMS.u, PARAMS.v];
    % symmetrization needed for directed case
    FlowMatrix = (FlowMatrix +FlowMatrix')/2;
    
else % undirected
    FlowMatrix = diag(PARAMS.pi)*PARAMS.P;
    PARAMS.null_vec = [PARAMS.pi, PARAMS.pi];
end

end
