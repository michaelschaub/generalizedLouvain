function [FlowMatrix, PARAMS] = louvain_FNL(Graph, time, PARAMS)
% Computes the full normalised stabilty matrix

PARAMS.louvain_type = 'normalised';
PARAMS.fixed_time = true;

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        dout = sum(Graph,2);
        dangling = (dout==0);
        dout(dangling) = 1;
        Dout = sparse(diag(dout));
        clear dout;
        M = (1-PARAMS.teleport_tau)*(Dout\Graph); % deterministic part of transition
        
        M =	M + diag(PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes)/PARAMS.NbNodes;
        
        clear Dout dangling
        
        if PARAMS.eig_decomp == false
            [v, lambda_all] = eigs(M'); % largest eigenvalue of transition matrix corresponds to stat.distribution.
            lambda = max(diag(real(lambda_all)));
            v = v(:,diag(lambda_all) == lambda);
            v = abs(v);              % make sure eigenvector is positive
            clear lambda;
            % store results for future use
            PARAMS.precomputed = true;
            PARAMS.pi = v/sum(v);
            PARAMS.P = M;
            
        elseif PARAMS.eig_decomp
            [V, lambda] = eig(full(M')); % M' = V*lambda/V
            lambda_max = max(diag(real(lambda)));
            v = V(:,diag(lambda) == lambda_max);
            v = abs(v);              % make sure eigenvector is positive
            % store results for future use
            PARAMS.precomputed = true;
            PARAMS.pi = v/sum(v);
            PARAMS.P = V;
            PARAMS.lambda = diag(lambda);
            
        else
            error('Something went wrong here')
        end
        
        if PARAMS.precompute_only
            FlowMatrix = 0;
            return
        end
        
        
        % undirected case
    else
        
        trans=sparse(diag( (sum(Graph)).^(-1) ) * Graph);
        % store results for future use
        PARAMS.precomputed = true;       
        PI=(diag(sum(Graph)))/sum(sum(Graph));  %diag matrix with stat distr
        PARAMS.pi = diag(PI);   % store results for future use
        
        if PARAMS.eig_decomp== false
            PARAMS.P = trans;
        elseif PARAMS.eig_decomp== true
            [V, lambda] = eig(full(trans')); % M' = V*lambda/V
            % normalise eigenvectors
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
    
end

% stationary distribution and "transition matrix" have been computed before
if PARAMS.directed == true
    if PARAMS.eig_decomp == true
        % We compute exp(M'-I) = V exp(lambda(M)-1) /V and the transpose
        % (P'(t)*PI)' == PI*P(t)
        FlowMatrix = ( (PARAMS.P*diag(exp((PARAMS.lambda-1)*time))...
            /PARAMS.P)*diag(PARAMS.pi) )';
    else
        FlowMatrix = diag(PARAMS.pi)*...
            expm(time* (PARAMS.P - eye(size(PARAMS.P))) );
    end
    % symmetrization needed for directed case
    FlowMatrix = real(FlowMatrix +FlowMatrix')/2;
    
else % undirected
    if PARAMS.eig_decomp == true
        % TODO: can use oblique orthogonality relations of M!?
        FlowMatrix = real( (PARAMS.P*diag(exp((PARAMS.lambda-1)*time))...
            /PARAMS.P)*diag(PARAMS.pi) );
    else
        FlowMatrix = diag(PARAMS.pi)*...
            expm(time* (PARAMS.P - eye(size(PARAMS.P))) );
    end
end

end
