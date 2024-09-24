function option_price = MonteCarloAmerican(S0, K, r, sigma, T, n, m, type)
    % S0 - Initial stock price
    % K - Strike price
    % r - Risk-free rate
    % sigma - Volatility
    % T - Time to maturity
    % n - Number of simulations
    % m - Number of time steps
    % type - Option type ('call' or 'put')
    
    % Time step
    dt = T / m;
    
    % Simulate asset paths
    S = zeros(n, m + 1);
    S(:, 1) = S0;
    for t = 2:m+1
        z = randn(n, 1);
        S(:, t) = S(:, t-1) .* exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * z);
    end
    
    % Calculate payoffs at maturity
    if type == "call"
        payoff = max(S(:, end) - K, 0);
    elseif type == "put"
        payoff = max(K - S(:, end), 0);
    end
    
    % Work backwards through the asset paths
    for t = m:-1:2
        % Indices of in-the-money paths
        if type == "call"
            itm = find(S(:, t) > K);
        elseif type == "put"
            itm = find(S(:, t) < K);
        end
        
        % Discounted payoffs from the next time step
        discounted_payoff = exp(-r * dt) * payoff(itm);
        
        % Regression on the in-the-money paths
        X = S(itm, t);
        Y = discounted_payoff;
        
        % Basis functions (polynomials up to degree 2)
        A = [ones(length(X), 1), X, X.^2];
        
        % Regression coefficients
        b = A \ Y;
        
        % Continuation values
        continuation_value = A * b;
        
        % Exercise values
        if type == "call"
            exercise_value = max(S(itm, t) - K, 0);
        elseif type == "put"
            exercise_value = max(K - S(itm, t), 0);
        end
        
        % Determine whether to exercise or continue
        exercise = exercise_value > continuation_value;
        
        % Update payoffs
        payoff(itm(exercise)) = exercise_value(exercise);
        payoff(itm(~exercise)) = discounted_payoff(~exercise);
    end
    
    % Discount the payoffs back to the present
    option_price = exp(-r * dt) * mean(payoff);
end
