%% On apple 
callsDataApple = readtable('apple_option_data.xlsx');
putsDataApple = readtable('AAPL_puts_option_data.xlsx');
callsDataApple.UnderlyingPrice = str2double(callsDataApple.UnderlyingPrice);
callsDataApple.OptionPrice = str2double(callsDataApple.OptionPrice);
putsDataApple.UnderlyingPrice = str2double(putsDataApple.UnderlyingPrice);
putsDataApple.OptionPrice = str2double(putsDataApple.OptionPrice);
% putsDataApple = putsDataApple((putsDataApple.strike<=205)&(putsDataApple.strike>=180), :);
% disp(putsDataApple);
callsDataApple.lastTradeDate = datetime(callsDataApple.lastTradeDate, 'InputFormat', 'yyyy-MM-dd');
putsDataApple.lastTradeDate = datetime(putsDataApple.lastTradeDate, 'InputFormat', 'yyyy-MM-dd');
% Remove the time zone to make them compatible
% putsDataApple.lastTradeDate.TimeZone = '';
% Assuming callsData is already loaded
r=0.05;
E=180;
% 
% % Initialize new arrays to hold the extracted numeric values
% extractedBetweenCommas = nan(height(callsDataApple), 1);
% extractedAfterLastComma = nan(height(callsDataApple), 1);
% 
% % Loop through each row to extract the values
% for i = 1:height(callsDataApple)
%     % Get the string from the second column
%     str = char(callsDataApple{i, 2});
%     
%     % Find the positions of the commas
%     commaIndices = strfind(str, ',');
%     
%     % Extract the value between the first two commas
%     if length(commaIndices) >= 2
%         % Extract the value between the first two commas
%         valueBetween = str(commaIndices(1)+1 : commaIndices(2)-1);
%         % Trim any leading or trailing whitespace
%         valueBetween = strtrim(valueBetween);
%         % Convert the extracted value to numeric
%         valueBetweenNumeric = str2double(valueBetween);
%     else
%         % If there are not enough commas, set the value to NaN
%         valueBetweenNumeric = NaN;
%     end
%     
%     % Store the extracted numeric value in the array
%     extractedBetweenCommas(i) = valueBetweenNumeric;
%     
%     % Extract the value after the last comma
%     if ~isempty(commaIndices)
%         % Extract the value after the last comma
%         valueAfterLast = str(commaIndices(end)+1:end);
%         % Trim any leading or trailing whitespace
%         valueAfterLast = strtrim(valueAfterLast);
%         % Convert the extracted value to numeric
%         valueAfterLastNumeric = str2double(valueAfterLast);
%     else
%         % If there is no comma, set the value to NaN
%         valueAfterLastNumeric = NaN;
%     end
%     
%     % Store the extracted numeric value in the array
%     extractedAfterLastComma(i) = valueAfterLastNumeric;
% end
% 
% % Add the new columns to the table
% callsDataApple.Underlying = extractedBetweenCommas;
% callsDataApple.OptionPrice = extractedAfterLastComma;
% 
% 
% % Display the updated table to verify the new columns
% disp(callsDataApple);
% % Get the unique values of the 'strike' column
% unique_strikes = unique(putsDataApple.strike);
% unique_dates = unique(putsDataApple.lastTradeDate);

% % Create a cell array to store subtables
% subtables = cell(length(unique_strikes), 1);
% returnsPuts=zeros(length(unique_strikes));
% sigma_annualPuts=zeros(length(unique_strikes),1);
% ExpDatePuts = zeros(height(putsDataApple),1)+datetime('2024-05-24','InputFormat', 'yyyy-MM-dd');
% TtMPuts=(1/252)*countBusinessDaysColumn(putsDataApple.lastTradeDate,ExpDatePuts);
% % Loop through each unique strike value and create a subtable
% for i = 1:length(unique_strikes)
%     % Get the current strike value
%     current_strike = unique_strikes(i);
%     
%     % Create a subtable with rows that have the current strike value
%     subtables{i} = putsDataApple(putsDataApple.strike == current_strike, :);
%     sigma_annualPuts(i) = (std(diff(log(subtables{i}.underlyingPrice))))*sqrt(252);
% end
% putsDataApple.Vol= repmat(sigma_annualPuts, height(unique_dates), 1);
% putsDataApple.TtM=TtMPuts;

%volatility?
% Compute daily returns
returnsCalls = diff(log(callsDataApple.UnderlyingPrice));
returnsPuts = diff(log(putsDataApple.UnderlyingPrice));
% Compute the standard deviation of the daily returns
sigma_dailyCalls = std(returnsCalls);
sigma_dailyPuts = std(returnsPuts);
% Annualize the standard deviation (assuming 252 trading days in a year)
sigma_annualCalls = sigma_dailyCalls * sqrt(252);
sigma_annualPuts= sigma_dailyPuts * sqrt(252);

ExpDate = zeros(height(callsDataApple),1)+datetime('2024-05-24','InputFormat', 'yyyy-MM-dd');
ExpDatePuts = zeros(height(putsDataApple),1)+datetime('2024-05-24','InputFormat', 'yyyy-MM-dd');

TtM=(1/252)*countBusinessDaysColumn(callsDataApple.lastTradeDate,ExpDate);
TtMputs=(1/252)*countBusinessDaysColumn(putsDataApple.lastTradeDate,ExpDatePuts);
% Initialize an array for theoretical prices
theoreticalPricesCK = zeros(height(callsDataApple), 1);
theoreticalPricesMC = zeros(height(callsDataApple), 1);

% Compute theoretical prices for each row
for i = 1:height(callsDataApple)
    S = callsDataApple.UnderlyingPrice(i);
    K = E;
    T_i = TtM(i);
    theoreticalPricesCK(i) = AmCK(S, K, r,T_i,sigma_annualCalls, 600,600,1.2,1e-6, 'call');
    theoreticalPricesMC(i) = MonteCarloAmerican(S, K, r,sigma_annualCalls,T_i, 10000,100, 'call');
    
end

% Add the theoretical prices to the table
callsDataApple.theoreticalPricesCK = theoreticalPricesCK;
callsDataApple.theoreticalPricesMC = theoreticalPricesMC;
% Plot the data
figure;
hold on;
%plot(callsDataApple.Var1, callsDataApple.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(callsDataApple.lastTradeDate, callsDataApple.OptionPrice, '-x', 'DisplayName', 'Option Price');
plot(callsDataApple.lastTradeDate, callsDataApple.theoreticalPricesCK, '-x', 'DisplayName', 'theoreticalPricesCK');
plot(callsDataApple.lastTradeDate, callsDataApple.theoreticalPricesMC, '-x', 'DisplayName', 'theoreticalPricesMC');
title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
hold off;

%PLOT DIFFERENCE
callsDataApple.absDiff = abs(callsDataApple.theoreticalPricesCK-callsDataApple.OptionPrice);
callsDataApple.relDiff=100.*callsDataApple.absDiff./callsDataApple.OptionPrice;
figure;
hold on;
plot(callsDataApple.lastTradeDate, callsDataApple.absDiff, '-x', 'DisplayName', 'absDiff');
plot(callsDataApple.lastTradeDate, callsDataApple.relDiff, '-x', 'DisplayName', 'relDiff');
title('Abs & Rel Differences');
xlabel('Date');
ylabel('Differences');
legend('show');
grid on;
hold off;
%plot(cal

%IV
%IV
IV_CK=zeros(height(callsDataApple),1);
IV_MC=zeros(height(callsDataApple),1);
%IMPLIED VOL. CK
for j=1:10
tic
for i = 1:height(callsDataApple)
    S = callsDataApple.UnderlyingPrice(i);
    T_i = TtM(i);
    V=callsDataApple.OptionPrice(i);
    %VMC=callsDataTTE.theoreticalPricesMC(i);
    errorFunctionCK = @(sigma) (AmCK(S, K, r,T_i,sigma, 800,800,1.2,1e-6, 'call') - V)^2;
    %errorFunctionMC = @(sigma) (MonteCarloAmerican(S, K, r,sigma,T_i, 100000,100, 'call') - V)^2;
    
    IV_CK(i)=fminbnd(errorFunctionCK, 0.05, 1);
    %IV_MC(i)=fminbnd(errorFunctionMC, 0.05, 1);
    
end
timpoIVMC=toc;
end
%timpoIVCK=toc;
IV_MCmean=mean([IV_MC(:,1) IV_MC(:,2) IV_MC(:,3) IV_MC(:,4) IV_MC(:,5) IV_MC(:,6) IV_MC(:,7) IV_MC(:,8) IV_MC(:,9) IV_MC(:,10)],2);
callsDataApple.IV_CK=IV_CK;
callsDataApple.IV_MC=IV_MC;

figure;
hold on;
% for j =1:10
%     plot(TtM, IV_MC(:,j), '-','Color',[0.5,0.5,0.5],'HandleVisibility', 'off');
% end
plot(TtM, callsDataApple.IV_CK, '-', 'DisplayName', 'IV-CK');
%plot(TtM, callsDataApple.IVOL_MID, '-', 'DisplayName', 'IVOL');
plot(TtM, callsDataApple.IV_MC, '-','DisplayName', 'IV-MCmean');

%plot(TtM, callsDataTTE.HistVol, 'DisplayName', 'Historical Volatulity of the unerlying');
title('Implied Volatility for K=180 Calls');
xlabel('Time to Maturity');
ylabel('Implied Volatility');
legend('show');
grid on;
ylim([0.05 0.6]);
hold off;


relDifIVck=abs(callsDataNOVO.IV_CK-callsDataNOVO.IVOL_MID)./callsDataNOVO.IVOL_MID;
relDifIVmc=abs(callsDataNOVO.IV_MC-callsDataNOVO.IVOL_MID)./callsDataNOVO.IVOL_MID;
figure;
hold on;
plot(TtM, relDifIVck, '-x', 'DisplayName', 'relDifIVck');
plot(TtM, relDifIVmc, '-x', 'DisplayName', 'relDifIVmc');
title('Relative Difference IV');
xlabel('Time to Maturity');
ylabel('Rel Diff');
legend('show');
grid on;
ylim([0 1]);
hold off;


%For puts:
% Initialize an array for theoretical prices
theoreticalPricesPutsCK = zeros(height(putsDataApple), 1);
theoreticalPricesPutsMC = zeros(height(putsDataApple), 1);

% Compute theoretical prices for each row
for i = 1:height(putsDataApple)
    S = putsDataApple.UnderlyingPrice(i);
    K = 180;
    T_i = TtMputs(i);
    theoreticalPricesPutsCK(i) = AmCK(S, K, r,T_i,sigma_annualPuts, 800,800,1.2,1e-6, 'put');
    theoreticalPricesPutsMC(i) = MonteCarloAmerican(S, K, r,sigma_annualPuts,T_i, 10000,100, 'put');
    
end

% Add the theoretical prices to the table
putsDataApple.theoreticalPricesPutsCK = theoreticalPricesPutsCK;
putsDataApple.theoreticalPricesPutsMC = theoreticalPricesPutsMC;
% Plot the data
figure;
hold on;
%plot(putsDataApple.Var1, callsDataApple.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(putsDataApple.lastTradeDate, putsDataApple.OptionPrice, '-x', 'DisplayName', 'Option Price');
plot(putsDataApple.lastTradeDate, putsDataApple.theoreticalPricesPutsCK, '-x', 'DisplayName', 'theoreticalPricesPutsCK');
plot(putsDataApple.lastTradeDate, putsDataApple.theoreticalPricesPutsMC, '-x', 'DisplayName', 'theoreticalPricesPutsMC');

title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
hold off;


%PLOT DIFFERENCE
putsDataApple.absDiff = abs(putsDataApple.theoreticalPricesPutsCK-putsDataApple.OptionPrice);
putsDataApple.relDiff=100.*putsDataApple.absDiff./putsDataApple.OptionPrice;
figure;
hold on;
plot(putsDataApple.lastTradeDate, putsDataApple.absDiff, '-x', 'DisplayName', 'absDiff');
plot(putsDataApple.lastTradeDate, putsDataApple.relDiff, '-x', 'DisplayName', 'relDiff');
title('Abs & Rel Differences');
xlabel('Date');
ylabel('Differences');
legend('show');
grid on;
hold off;


IV_CK=zeros(height(putsDataApple),1);
IV_MC=zeros(height(putsDataApple),1);
%IMPLIED VOL. CK
for j=1:10
tic
for i = 1:height(putsDataApple)
    S = putsDataApple.UnderlyingPrice(i);
    T_i = TtMputs(i);
    V=putsDataApple.OptionPrice(i);
    %VMC=callsDataTTE.theoreticalPricesMC(i);
    %errorFunctionCK = @(sigma) (AmCK(S, K, r,T_i,sigma, 800,800,1.2,1e-6, 'put') - V)^2;
    errorFunctionMC = @(sigma) (MonteCarloAmerican(S, K, r,sigma,T_i, 100000,100, 'put') - V)^2;
    
    %IV_CK(i)=fminbnd(errorFunctionCK, 0.05, 1);
    IV_MC(i)=fminbnd(errorFunctionMC, 0.05, 1);
    
end
timpoIVMC=toc;
end
%timpoIVCK=toc;
IV_MCmean=mean([IV_MC(:,1) IV_MC(:,2) IV_MC(:,3) IV_MC(:,4) IV_MC(:,5) IV_MC(:,6) IV_MC(:,7) IV_MC(:,8) IV_MC(:,9) IV_MC(:,10)],2);
putsDataApple.IV_CK=IV_CK;
putsDataApple.IV_MC=IV_MC;

figure;
hold on;
plot(TtM, callsDataApple.IV_CK, '-', 'DisplayName', 'IV-CK');
plot(TtMputs, putsDataApple.IV_CK, '-', 'DisplayName', 'IV-CKputs');
%plot(TtM, callsDataApple.IVOL_MID, '-', 'DisplayName', 'IVOL');
plot(TtM, callsDataApple.IV_MC, '-','DisplayName', 'IV-MC');
plot(TtMputs, putsDataApple.IV_MC, '-','DisplayName', 'IV-MCputs');
%plot(TtM, callsDataTTE.HistVol, 'DisplayName', 'Historical Volatulity of the unerlying');
title('Implied Volatility for K=180 ');
xlabel('Time to Maturity');
ylabel('Implied Volatility');
legend('show');
grid on;
ylim([0.05 1]);
hold off;


%% IV Apple 31052024:
callsApple = readtable('AAPL_calls_option_data310524.xlsx');

callsApple.Strike=str2double(callsApple.Strike);
callsApple.OptionPrice=str2double(callsApple.OptionPrice);
callsApple.impliedVolatility=str2double(callsApple.impliedVolatility);
callsApple.UnderlyingPrice=str2double(callsApple.UnderlyingPrice);
callsApple=callsApple(callsApple.Strike>=150,:);
ExpDate = zeros(height(callsApple),1)+datetime('2024-05-31','InputFormat', 'yyyy-MM-dd');
TtM=(1/252)*countBusinessDaysColumn(callsApple.lastTradeDate,ExpDate);
callsApple.TimeToMaturity=TtM;
r=0.05;

uniqueStrikes=unique(callsApple.Strike);
ssigma=zeros(height(uniqueStrikes),1);
for i=1:height(uniqueStrikes)
    K=uniqueStrikes(i);
    underlyings=callsApple(callsApple.Strike==K,:).UnderlyingPrice;
    ssigma(i)= std(diff(log(underlyings)))* sqrt(252);
end
callsApple.sigma=zeros(height(callsApple),1);
for i=1:height(uniqueStrikes)
    K=uniqueStrikes(i);
    callsApple(callsApple.Strike==K,:).sigma=ssigma(i)+zeros(height(callsApple(callsApple.Strike==K,:).sigma),1);
    
end
%SSigma=repmat(ssigma, height(callsApple)/height(uniqueStrikes), 1);
theoreticalPricesCK = zeros(height(callsApple), 1);
theoreticalPricesMC = zeros(height(callsApple), 1);
estimated_IV_CK=zeros(height(callsApple), 1);
estimated_IV_MC=zeros(height(callsApple), 1);
for i=1:height(callsApple)
    S=callsApple.UnderlyingPrice(i);
    K=callsApple.Strike(i);
    T=callsApple.TimeToMaturity(i);
    sigma=callsApple.sigma(i);
    %theoreticalPricesCK(i)=AmCK(S,K,r,T,sigma,800,1000,1.2,1e-6,"call");
    %theoreticalPricesMC(i)=MonteCarloAmerican(S, K, r,sigma,T, 10000,100, 'call');
    
    V=callsApple.OptionPrice(i);
    %errorFunctionCK = @(x) (AmCK(S,K,r,T,x,600,600,1.2,1e-6,"call") - V)^2;
    errorFunctionMC = @(x) (MonteCarloAmerican(S, K, r,x,T, 100000,100, 'call') - V)^2;
    %estimated_IV_CK(i)=fminbnd(errorFunctionCK, 0.05, 1.5);
    estimated_IV_MC(i)=fminbnd(errorFunctionMC, 0.05, 1.5);
end
callsApple.estimated_IV_CK=estimated_IV_CK;
callsApple.estimated_IV_MC=estimated_IV_MC;
callsApple.theoreticalPricesMC=theoreticalPricesMC;
callsApple.theoreticalPricesCK=theoreticalPricesCK;

[X, Y] = meshgrid(uniqueStrikes, unique(TtM));

% Initialize Z (IV values) with NaNs

Z=table2array(estimatedIVMC_Elisa)';
Z=estimatedIVMC_Elisa';
ZCN=zeros(height(unique(TtM)),height(uniqueStrikes));
ZMC=zeros(height(unique(TtM)),height(uniqueStrikes));
for i=1:height(uniqueStrikes)
    K=uniqueStrikes(i);
    ZCN(:,i)=callsApple(callsApple.Strike==K,:).estimated_IV_CK;
    ZMC(:,i)=callsApple(callsApple.Strike==K,:).estimated_IV_MC;
end

% Create the 3D surface plot
figure;
surf(X, Y, ZCN);

% Plot the data
figure;
hold on;
%plot(putsDataApple.Var1, callsDataApple.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(callsApple.lastTradeDate, callsApple.impliedVolatility, '-x', 'DisplayName', 'Option Price');
plot(callsApple.lastTradeDate, callsApple.estimated_IV_CK, '-x', 'DisplayName', 'theoreticalPricesCK');
plot(callsApple.lastTradeDate, callsApple.estimated_IV_MC, '-x', 'DisplayName', 'theoreticalPricesMC');

title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
hold off;

%plot one line of time to maturity:
figure;
hold on;
% Plot Real IV for each column
for i=7:height(uniqueStrikes)
    K=uniqueStrikes(i);
    plot(callsApple(callsApple.Strike==K,:).TimeToMaturity, callsApple(callsApple.Strike==K,:).impliedVolatility, '-', 'DisplayName', ['Real IV T= ']);
    plot(callsApple(callsApple.Strike==K,:).TimeToMaturity, callsApple(callsApple.Strike==K,:).estimated_IV_CK, '-x', 'DisplayName', ['estimated_IV_CK ']);
    plot(callsApple(callsApple.Strike==K,:).TimeToMaturity, callsApple(callsApple.Strike==K,:).estimated_IV_MC, '-x', 'DisplayName', ['estimated_IV_MC ']);
end

title('Real IV Data vs Computed Crank Nicolson');
xlabel('Strike');
ylabel('IV');
% legend('show');
grid on;
hold off

%3D
%Impliet Volatility 3D surface
% Extract relevant columns
% Create a grid for strike and time to maturity
[X, Y] = meshgrid(uniqueStrikes, unique(TtM));

% Initialize Z (IV values) with NaNs
Zck = zeros(size(X));

% Fill Z with the IV values corresponding to the grid points
for i = 1:height(callsApple)
    % Find the indices in the grid
    xi = find(X(1,:) == callsApple.Strike(i));
    yi = find(Y(:,1) == callsApple.TimeToMaturity(i));
    
    % Assign the IV value
    Zck(yi, xi) = callsApple.estimated_IV_CK(i);
end
Z=flipud(Z);
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

% Replace 0 values with the value to the left
[m, n] = size(Z); % Get the size of the matrix

for i = 1:m
    for j = 2:n % Start from the second column
        if Z(i, j) == 0
            Z(i, j) = Z(i, j-1);
        end
    end
end
writetable(callsApple,'callsApple.txt')