function [FlowMatrix, PARAMS] = louvain_FCL(Graph, time, PARAMS)
% Computes the full combinatorial stabilty matrix

PARAMS.louvain_type = 'normalised';
PARAMS.fixed_time = true;

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        dout = sum(Graph,2);
        dangling = (dout==0);
        dout(dangling) = 1;
        av_degree = mean(dout);
        Dout = sparse(diag(dout));
        clear dout;
        M = (1-PARAMS.teleport_tau)*(Dout\Graph); % deterministic part of transition
        
        M =	M + diag(PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes)/PARAMS.NbNodes;
        L= (Dout./av_degree)*(eye(size(M))- M);
        clear Dout dangling
        
        if PARAMS.eig_decomp == false
            % get eigenvalues and eigenvectors of Laplacian
            [v, lambda_all] = eigs(L',5,'SR'); 
            lambda = min(diag(real(lambda_all)));

            % get stationary dist and normalize it
            v = abs(v(:,diag(lambda_all) == lambda));
            v = v/sum(v);             
            clear lambda;
            
            % store results for future use
            PARAMS.precomputed = true;
            PARAMS.pi = v;
            PARAMS.P = -L;
            
        elseif PARAMS.eig_decomp
            [V, lambda] = eig(full(L')); % L' = V*lambda/V
            lambda0 = min(diag(real(lambda)));
            v = V(:,diag(lambda) == lambda0);
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
        degree = sum(Graph,2);
        av_degree = mean(degree);
        L = (diag(degree)-Graph)/av_degree;
        % store results for future use
        PARAMS.precomputed = true;       
        PI=eye(PARAMS.NbNodes)/PARAMS.NbNodes;  %diag matrix with stat distr
        PARAMS.pi = diag(PI);   % store results for future use
        
        if PARAMS.eig_decomp== false
            PARAMS.P = -L;
        elseif PARAMS.eig_decomp== true
            [V, lambda] = eig(full(L)); % M' = V*lambda/V
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
        % L'=V*Lambda/V we are computing exp(L') and then transpose!
        % (P'(t)*PI)' == PI*P(t)
        FlowMatrix = ( (PARAMS.P*diag(exp((-PARAMS.lambda)*time))...
            /PARAMS.P)*diag(PARAMS.pi) )';
    else
        FlowMatrix = diag(PARAMS.pi)*...
            expm(time* (PARAMS.P));
    end
    % symmetrization needed for directed case
    FlowMatrix = real(FlowMatrix +FlowMatrix')/2;
    
else % undirected
    if PARAMS.eig_decomp == true
        % Same as above but now we know that P'(t)*PI == PI*P(t), also we know that L is symmetric,
        % so we can use transpose of eigenvectors
        FlowMatrix = real( (PARAMS.P*diag(exp((-PARAMS.lambda)*time))...
            *PARAMS.P')*diag(PARAMS.pi) );
    else
        FlowMatrix = diag(PARAMS.pi)*...
            expm(time* (PARAMS.P)) ;
    end
end

end
