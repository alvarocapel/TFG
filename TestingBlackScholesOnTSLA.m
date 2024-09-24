%% On TSlA
callsDataTSLA = readtable('tsla_calls_option_data.xlsx');
callsDataTSLA.UnderlyingPrice = str2double(callsDataTSLA.UnderlyingPrice);
callsDataTSLA.OptionPrice = str2double(callsDataTSLA.OptionPrice);
% putsDataTSLA = putsDataTSLA((putsDataTSLA.strike<=205)&(putsDataTSLA.strike>=180), :);
% disp(putsDataTSLA);
callsDataTSLA.lastTradeDate = datetime(callsDataTSLA.lastTradeDate, 'InputFormat', 'yyyy-MM-dd');
% Remove the time zone to make them compatible
% putsDataTSLA.lastTradeDate.TimeZone = '';
% Assuming callsData is already loaded
r=0.05;
E=170;

%volatility?
% Compute daily returns
returnsCalls = diff(log(callsDataTSLA.UnderlyingPrice));
% Compute the standard deviation of the daily returns
sigma_dailyCalls = std(returnsCalls);
% Annualize the standard deviation (assuming 252 trading days in a year)
sigma_annualCalls = sigma_dailyCalls * sqrt(252);


ExpDate = zeros(height(callsDataTSLA),1)+datetime('2024-05-31','InputFormat', 'yyyy-MM-dd');

TtM=(1/252)*countBusinessDaysColumn(callsDataTSLA.lastTradeDate,ExpDate);
% Initialize an array for theoretical prices
theoreticalPricesCK = zeros(height(callsDataTSLA), 1);
theoreticalPricesMC = zeros(height(callsDataTSLA), 1);

% Compute theoretical prices for each row
for i = 1:height(callsDataTSLA)
    S = callsDataTSLA.UnderlyingPrice(i);
    K = E;
    T_i = TtM(i);
    theoreticalPricesCK(i) = AmCK(S, K, r,T_i,sigma_annualCalls, 600,600,1.2,1e-6, 'call');
    theoreticalPricesMC(i) = MonteCarloAmerican(S, K, r,sigma_annualCalls,T_i, 100000,100, 'call');
    
end

% Add the theoretical prices to the table
callsDataTSLA.theoreticalPricesCK = theoreticalPricesCK;
callsDataTSLA.theoreticalPricesMC = theoreticalPricesMC;
% Plot the data
figure;
hold on;
%plot(callsDataTSLA.Var1, callsDataTSLA.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(callsDataTSLA.lastTradeDate, callsDataTSLA.OptionPrice, '-x', 'DisplayName', 'Option Price');
plot(callsDataTSLA.lastTradeDate, callsDataTSLA.theoreticalPricesCK, '-x', 'DisplayName', 'theoreticalPricesCK');
plot(callsDataTSLA.lastTradeDate, callsDataTSLA.theoreticalPricesMC, '-x', 'DisplayName', 'theoreticalPricesMC');
title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
ylim([0 50]);
hold off;

% %PLOT DIFFERENCE
% callsDataTSLA.absDiff = abs(callsDataTSLA.theoreticalPricesCK-callsDataTSLA.OptionPrice);
% callsDataTSLA.relDiff=100.*callsDataTSLA.absDiff./callsDataTSLA.OptionPrice;
% figure;
% hold on;
% plot(callsDataTSLA.lastTradeDate, callsDataTSLA.absDiff, '-x', 'DisplayName', 'absDiff');
% plot(callsDataTSLA.lastTradeDate, callsDataTSLA.relDiff, '-x', 'DisplayName', 'relDiff');
% title('Abs & Rel Differences');
% xlabel('Date');
% ylabel('Differences');
% legend('show');
% grid on;
% hold off;

%Implied volatilty:




%% TSLA puts
putsDataTSLA = readtable('tsla_puts_option_data.xlsx');

putsDataTSLA.UnderlyingPrice = str2double(putsDataTSLA.UnderlyingPrice);
putsDataTSLA.OptionPrice = str2double(putsDataTSLA.OptionPrice);

putsDataTSLA.lastTradeDate = datetime(putsDataTSLA.lastTradeDate, 'InputFormat', 'yyyy-MM-dd');
r=0.05;
E=170;
%volatility?
% Compute daily returns

returnsPuts = diff(log(putsDataTSLA.UnderlyingPrice));
% Compute the standard deviation of the daily returns

sigma_dailyPuts = std(returnsPuts);
% Annualize the standard deviation (assuming 252 trading days in a year)

sigma_annualPuts= sigma_dailyPuts * sqrt(252);
ExpDatePuts = zeros(height(putsDataTSLA),1)+datetime('2024-05-31','InputFormat', 'yyyy-MM-dd');
TtMputs=(1/252)*countBusinessDaysColumn(putsDataTSLA.lastTradeDate,ExpDatePuts);

%For puts:
% Initialize an array for theoretical prices
theoreticalPricesPutsCK = zeros(height(putsDataTSLA), 1);
theoreticalPricesPutsMC = zeros(height(putsDataTSLA), 1);

% Compute theoretical prices for each row
for i = 1:height(putsDataTSLA)
    S = putsDataTSLA.UnderlyingPrice(i);
    K = E;
    T_i = TtMputs(i);
    theoreticalPricesPutsCK(i) = AmCK(S, K, r,T_i,sigma_annualPuts, 600,600,1.2,1e-6, 'put');
    theoreticalPricesPutsMC(i) = MonteCarloAmerican(S, K, r,sigma_annualPuts,T_i, 100000,100, 'put');
    
end
% Add the theoretical prices to the table
putsDataTSLA.theoreticalPricesPutsCK = theoreticalPricesPutsCK;
putsDataTSLA.theoreticalPricesPutsMC = theoreticalPricesPutsMC;
% Plot the data
figure;
hold on;
%plot(putsDataTSLA.Var1, callsDataTSLA.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(putsDataTSLA.lastTradeDate, putsDataTSLA.OptionPrice, '-x', 'DisplayName', 'Option Price');
plot(putsDataTSLA.lastTradeDate, putsDataTSLA.theoreticalPricesPutsCK, '-x', 'DisplayName', 'theoreticalPricesPutsCK');
plot(putsDataTSLA.lastTradeDate, putsDataTSLA.theoreticalPricesPutsMC, '-x', 'DisplayName', 'theoreticalPricesPutsMC');

title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
ylim([0 45]);
hold off;


% %PLOT DIFFERENCE
% putsDataTSLA.absDiff = abs(putsDataTSLA.theoreticalPricesPutsCK-putsDataTSLA.OptionPrice);
% putsDataTSLA.relDiff=100.*putsDataTSLA.absDiff./putsDataTSLA.OptionPrice;
% figure;
% hold on;
% plot(putsDataTSLA.lastTradeDate, putsDataTSLA.absDiff, '-x', 'DisplayName', 'absDiff');
% plot(putsDataTSLA.lastTradeDate, putsDataTSLA.relDiff, '-x', 'DisplayName', 'relDiff');
% title('Abs & Rel Differences');
% xlabel('Date');
% ylabel('Differences');
% legend('show');
% grid on;
% hold off;


%Implied volatilty:
putsDataTSLAIV = readtable('tsla_puts_option_dataIV.xlsx');

putsDataTSLAIV.UnderlyingPrice = str2double(putsDataTSLAIV.UnderlyingPrice);
putsDataTSLAIV.OptionPrice = str2double(putsDataTSLAIV.OptionPrice);
putsDataTSLAIV.Strike = str2double(putsDataTSLAIV.Strike);
putsDataTSLAIV.lastTradeDate = datetime(putsDataTSLAIV.lastTradeDate, 'InputFormat', 'yyyy-MM-dd');
r=0.05;

ExpDatePuts = zeros(height(putsDataTSLAIV),1)+datetime('2024-05-31','InputFormat', 'yyyy-MM-dd');
TtMputs=(1/252)*countBusinessDaysColumn(putsDataTSLAIV.lastTradeDate,ExpDatePuts);
putsDataTSLAIV.TimeToMaturity=TtMputs;
% Error function to minimize
for i = 1:height(putsDataTSLAIV)
    V=putsDataTSLAIV.OptionPrice(i);
    S=putsDataTSLAIV.UnderlyingPrice(i);
    E=putsDataTSLAIV.Strike(i);
    TimeToMaturity=putsDataTSLAIV.TimeToMaturity(i);
    %errorFunction = @(sigma) (AmCK(S, E, r,TimeToMaturity, sigma,800,800,1.2,1e-6, 'put') - V)^2;
    errorFunction = @(sigma) (MonteCarloAmerican(S, E, r, sigma,TimeToMaturity,100000,100,'put') - V)^2;
    
    % Use Brent's method to find the implied volatility
    % fminbnd is the MATLAB equivalent to optimize.brent in Python
    estimatedIV = fminbnd(errorFunction, 0.05, 1);
    putsDataTSLAIV.IV_MC(i)=estimatedIV;
end

%Impliet Volatility 3D surface
% Extract relevant columns
% Create a grid for strike and time to maturity
[X, Y] = meshgrid(unique(putsDataTSLAIV.Strike), unique(putsDataTSLAIV.TimeToMaturity));

% Initialize Z (IV values) with NaNs
Z = zeros(size(X));

% Fill Z with the IV values corresponding to the grid points
for i = 1:length(putsDataTSLAIV.Strike)
    % Find the indices in the grid
    xi = find(X(1,:) == putsDataTSLAIV.Strike(i));
    yi = find(Y(:,1) == putsDataTSLAIV.TimeToMaturity(i));
    
    % Assign the IV value
    Z(yi, xi) = putsDataTSLAIV.IV_MC(i);
end

figure;
% Number of subplots based on the columns in Z
numCols = size(Z, 2);
% Define the titles for each subplot
titles = {
    'Implied Volatility for K=160 Calls',
    'Implied Volatility for K=165 Calls',
    'Implied Volatility for K=170 Calls',
    'Implied Volatility for K=175 Calls',
    'Implied Volatility for K=180 Calls'
};
% Plot each column of Z in a separate subplot
for col = 1:numCols
    if col <= 3
        subplot(3, 2, col); % Top row, 1st to 3rd position in a 3x2 grid
    else
        subplot(3, 2, col + 1); % Bottom row, 4th and 5th position in a 3x2 grid
    end
    plot(Y(:,1), Z(:,col), '-', 'DisplayName', ['IV-CK for K=' num2str(160 + 5 * (col - 1))]);
    title(titles{col});
    xlabel('Time to Maturity');
    ylabel('Implied Volatility');
    grid on;
    legend show;
end

% Adjust layout
sgtitle('Implied Volatility for Different Strikes'); % Super title for the entire figure

% Create the 3D surface plot
figure;
surf(X, Y, Z);

% Add labels and title
xlabel('Strike Price');
ylabel('Time to Maturity');
zlabel('Implied Volatility');
title('Implied Volatility Surface');

% Improve visualization
colorbar;
shading interp;
% Extract the first slice of Z where TimeToMaturity is at its first value
firstSlice = Z(1, :);

% Create a 2D plot for this slice
figure;
plot(X(1, :), firstSlice, '-o');

% Add labels and title
xlabel('Strike Price');
ylabel('Implied Volatility');
title(['Implied Volatility at Time to Maturity = ', num2str(Y(1, 1))]);

grid on;
