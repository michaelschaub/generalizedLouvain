function [FlowMatrix, PARAMS] = louvain_modularity(Graph, time, PARAMS)


if PARAMS.directed
    PARAMS.louvain_type = 'generalised';
else
    PARAMS.louvain_type = 'normalised';
end
PARAMS.fixed_time = true;

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        dout = sum(Graph,2);
        din = sum(Graph,1);
        tot_w = sum(dout);
        PARAMS.null_vec = [dout/tot_w, din'/tot_w];
        PARAMS.P = Graph/tot_w;
        if PARAMS.precompute_only
            FlowMatrix = 0;
            return
        end
        
        
        
    else % undirected case
       
        
        if PARAMS.precompute_only
            FlowMatrix = 0;
            return
        end
        
    end
    
end

% stationary distribution and "transition matrix" have been computed before
if PARAMS.directed == true
    FlowMatrix = PARAMS.P;
    % symmetrization needed for directed case, null model takes care of
    % dir aspects
    FlowMatrix = (FlowMatrix +FlowMatrix')/2;
    
else % undirected
    FlowMatrix = Graph;
end

end
