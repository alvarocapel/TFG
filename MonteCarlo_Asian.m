function option_price = MonteCarlo_Asian(S0, K, r, sigma, T, n)
    % Generate a matrix of normal nx101 random numbers
    z = randn(n, 100);

    % Create matrices to store Brownian motion and asset prices
    brownian = zeros(n, 101);
    S = zeros(n, 101);

    % Create a vector to store the average asset prices over the Asian option period
    avg_S = zeros(n, 1);

    % Create a vector to store the payoff
    payoff = zeros(n, 1);

    % Set the initial value of the stock to the first column of the Asset Prices matrix
    S(:, 1) = S0;

    % For each simulation
    for i = 1:n
        % For each timestep (we split T in 101 times)
        for j = 2:101
            % Generate a Brownian motion of the timestep j*T/101
            brownian(i, j) = brownian(i, j - 1) + sqrt(T / 101) * z(i, j - 1);
            
            % Generate the asset path
            S(i, j) = S0 * exp((r - 0.5 * sigma^2) * (j * T / 101) + sigma * brownian(i, j));
        end

        % Calculate the average asset price over the Asian option period
        avg_S(i) = mean(S(i, 2:end));

        % Calculate the payoff of the Asian option
        payoff(i) = max(avg_S(i) - K, 0);
    end

    % Discount the payoff to move from T to now
    option_price = exp(-r * T) * mean(payoff);
end
