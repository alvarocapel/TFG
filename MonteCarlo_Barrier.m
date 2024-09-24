function option_price = MonteCarlo_Barrier(S0, K, r, sigma, T, b, n)
    % Generate a matrix of normal nx101 random numbers
    z = randn(n, 101); %we divide the time into 100 steps

    % Create matrices to store Brownian motion and asset prices
    brownian = zeros(n, 101);
    S = zeros(n, 101);

    % Create a vector to store the payoff
    payoff = zeros(n, 1);

    % Set the initial value of the stock to the first column of the Asset Prices matrix
    S(:, 1) = S0;
    countOfReached=0;
    % For each simulation
    for i = 1:n
        % For each timestep (we split T in 101 times)
        for j = 2:101
            % Generate a Brownian motion of the timestep j*T/101
            brownian(i, j) = brownian(i, j - 1) + sqrt(T / 101) * z(i, j - 1);
            
            % Generate the asset path
            S(i, j) = S0 * exp((r - 0.5 * sigma^2) * (j * T / 101) + sigma * brownian(i, j));
        end

        % Once we have the full path, check if the barrier is reached
        % if so, the option is deactivated
        if max(S(i, :)) > b
            payoff(i) = 0;
            countOfReached=countOfReached+1;
            plott=i;
        % If not, we have a European option
        else
            payoff(i) = max(S(i, end) - K, 0);
        end
    end

    % Discount the payoff to move from T to now
    option_price = exp(-r * T) * mean(payoff);
    plot(S(plott,:))
end
