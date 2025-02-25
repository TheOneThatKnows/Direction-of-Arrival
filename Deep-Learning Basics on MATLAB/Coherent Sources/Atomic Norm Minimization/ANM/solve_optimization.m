function [W, u, R0] = solve_optimization(G, R_hat)
    % Input:
    % G - matrix
    % R_hat - matrix
    % Output:
    % W - optimized matrix
    % u - vector for Toeplitz matrix
    % R0 - optimized matrix
    
    [N, M] = size(G); % Get dimension from G matrix
    mu = 1; % selected as 1 in the article
    
    cvx_begin sdp
        % Define variables
        variable W(M,M) symmetric
        variable u(N,1)
        variable R0(N,M)

        % Objective function
        objective = (mu/(2*sqrt(N)))*(trace(W) + trace(toeplitz(u'))) + ...
            0.5 * sum(sum((R0 .* G - R_hat) .* conj(R0 .* G - R_hat)))
            % 0.5*norm(R0 .* G - R_hat, 'fro')^2;
        
        minimize(objective)
        
        % Constraints
        % Positive semidefinite constraint
        [W, R0'; R0, toeplitz(u')] >= 0;
        
    cvx_end
end