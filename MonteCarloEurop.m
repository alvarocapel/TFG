%% Monte Carlo methods:
function option_price = MonteCarloEurop(S0, K, r, sigma, T, n,type)
    % Initialize vectors
    payoff = zeros(n, 1);
    S = zeros(n, 1);

    % Generate n different random numbers following a N(0,1)
    z = randn(n, 1);

    % Calculate the stock value at maturity (it is a random variable)
    S = S0 * exp(T * (r - sigma^2 / 2) + sigma * sqrt(T) * z);

    % Calculate the payoff
    if type=="call"
        payoff = max(S - K, 0);
    elseif type== "put" 
        payoff = max(K-S, 0);
    end
    % Return the discounted payoff
    option_price = exp(-r * T) * mean(payoff);
end
