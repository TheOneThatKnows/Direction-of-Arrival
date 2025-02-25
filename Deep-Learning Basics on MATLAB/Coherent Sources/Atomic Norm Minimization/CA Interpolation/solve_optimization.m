% does not work
function F_tilda = solve_optimization(C, F_hat)
    % Input:
    % C - binary matrix
    % F_hat - covariance matrix with zero diagonals
    % Output:
    % F_tilda - optimized covariance matrix
    
    N2 = size(F_hat, 1); % Get dimension from G matrix
    
    e = rand;
    cvx_begin
        % Define variables
        variable F_tilda(N2,N2)

        % Objective function
        objective = norm_nuc(F_tilda)
        
        minimize(objective)
        
        % Constraints
        subject to
            norm(F_tilda .* C - F_hat, 'fro') <= e
    cvx_end
end