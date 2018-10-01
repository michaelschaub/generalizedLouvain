function [FlowMatrix, PARAMS] = louvain_LNL(Graph, time, PARAMS)
% Computes the linearised normalised stabilty matrix

if PARAMS.directed
    PARAMS.louvain_type = 'generalised';
else
    PARAMS.louvain_type = 'normalised';
end
PARAMS.fixed_time = false;

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        dout = sum(Graph,2);
        dangling = (dout==0);
        dout(dangling) = 1;
        
        Dout = sparse(1:PARAMS.NbNodes,1:PARAMS.NbNodes,dout);
        clear dout;
        M1 = (1-PARAMS.teleport_tau)*(Dout\Graph); % deterministic part of transition
        
        a = sparse(1:PARAMS.NbNodes,1:PARAMS.NbNodes,PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes,1);
        b = ones(PARAMS.NbNodes,1)/PARAMS.NbNodes;

        
        % M = M1 + a*b';
        clear Dout dangling
        
        
        f = @(x) M1'*x + b*(a'*x);
        [V , lambda]=eigs(f,PARAMS.NbNodes);
        % BEWARE, this might be slow!
        %[V, lambda] = eigs(full(M')); % M' = V*lambda/V
        lambda_max = max(diag(real(lambda)));
        v = V(:,diag(lambda) == lambda_max);
        v = abs(v);              % make sure eigenvector is positive
        % store results for future use
        PARAMS.precomputed = true;
        PARAMS.pi = v/sum(v);
        PI = sparse(1:PARAMS.NbNodes,1:PARAMS.NbNodes,PARAMS.pi);
        PARAMS.P = PI*M1;
        
        PARAMS.u = -PI*a;
        PARAMS.v = b;
        
        if PARAMS.precompute_only
            FlowMatrix = 0;
            return
        end
        
        
        
    else % undirected case
        
        trans=sparse(diag( (sum(Graph)).^(-1) ) * Graph);
        % store results for future use
        PARAMS.precomputed = true;       
        PI=(diag(sum(Graph)))/sum(sum(Graph));  %diag matrix with stat distr
        PARAMS.pi = diag(PI);   % store results for future use      
        PARAMS.P = trans;    
        
        if PARAMS.precompute_only
            FlowMatrix = 0;
            return
        end
        
    end
    
end

% stationary distribution and "transition matrix" have been computed before
if PARAMS.directed == true
    FlowMatrix = PARAMS.P;
    PARAMS.null_vec = [PARAMS.pi, PARAMS.pi, time*PARAMS.u, PARAMS.v];
    % symmetrization needed for directed case
    FlowMatrix = (FlowMatrix +FlowMatrix')/2;
    
else % undirected
    FlowMatrix = sparse(1:PARAMS.NbNodes,1:PARAMS.NbNodes,PARAMS.pi)*PARAMS.P;
end

end
