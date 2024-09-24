%% ON Airbus 
%puts
callsDataAIR = readtable('Dades Opcions.xlsx','Sheet','AIRBUS');
disp(callsDataAIR);
rowsWithNaN = any(ismissing(callsDataAIR), 2);
% Remove rows with NaN values
callsDataAIR = callsDataAIR(~rowsWithNaN, :);
callsDataAIRIVputs = readtable('Dades Opcions.xlsx','Sheet','AIRBUS').IVOL_MIDputs;
rowsWithNaN = any(ismissing(callsDataAIRIVputs), 2);
% Remove rows with NaN values
callsDataAIRIVputs = callsDataAIRIVputs(~rowsWithNaN, :);

callsDataAIR.IVOL_MIDputs=callsDataAIRIVputs;
callsDataAIR.IVOL_MIDputs=(1/100)*callsDataAIR.IVOL_MIDputs;
K=68;
ExpDate = datetime('2024-06-21','InputFormat', 'yyyy-MM-dd');
r=0.025; %German risk free rate
%volatility?
% Compute daily returns
returnsPuts = diff(log(callsDataAIR.Underlying));
% Compute the standard deviation of the daily returns
sigma_dailyPuts = std(returnsPuts);
% Annualize the standard deviation (assuming 252 trading days in a year)
sigma_annualPuts = sigma_dailyPuts * sqrt(252);

theoreticalPricesCK=zeros(height(callsDataAIR),1);
theoreticalPricesMC=zeros(height(callsDataAIR),1);
% Compute theoretical prices for each row
for i = 1:height(callsDataAIR)
    S = callsDataAIR.Underlying(i);
    T_i =callsDataAIR.TimeToMaturity(i);
    K=callsDataAIR.StrikPuts(i);
    theoreticalPricesCK(i) = CK_Europ(S, K, r,sigma_annualPuts,T_i, 800,1500, 'put');
    theoreticalPricesMC(i) = MonteCarloEurop(S, K, r,sigma_annualPuts,T_i, 10000, 'put');
    
end

% Add the theoretical prices to the table
callsDataAIR.theoreticalPricesCK = theoreticalPricesCK;
callsDataAIR.theoreticalPricesMC = theoreticalPricesMC;

% Plot the data
figure;
hold on;
%plot(putsDataApple.Var1, callsDataApple.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(callsDataAIR.Date, callsDataAIR.Puts, '-x', 'DisplayName', 'Option Price');
plot(callsDataAIR.Date, callsDataAIR.theoreticalPricesCK, '-x', 'DisplayName', 'theoreticalPricesCK');
plot(callsDataAIR.Date, callsDataAIR.theoreticalPricesMC, '-x', 'DisplayName', 'theoreticalPricesMC');

title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
hold off;

IV_CK=zeros(height(callsDataAIR),1);
IV_MC=zeros(height(callsDataAIR),10);

%IMPLIED VOL. CK
for j=1:10
for i = 1:height(callsDataAIR)
    S = callsDataAIR.Underlying(i);
    T_i = callsDataAIR.TimeToMaturity(i);
    V=callsDataAIR.Puts(i);
    K=callsDataAIR.StrikPuts(i);
    %VMC=callsDataTTE.theoreticalPricesMC(i);
    %errorFunctionCK = @(sigma) (CK_Europ(S, K, r, sigma,T_i,800,1500, 'put') - V)^2;
    errorFunctionMC = @(sigma) (MonteCarloEurop(S, K, r, sigma,T_i,100000, 'put') - V)^2;
    
    %IV_CK(i)=fminbnd(errorFunctionCK, 0.1, 1);
    IV_MC(i,j)=fminbnd(errorFunctionMC, 0.1, 1);
    
end
end
IV_MCmean=mean([IV_MC(:,1) IV_MC(:,2) IV_MC(:,3) IV_MC(:,4) IV_MC(:,5) IV_MC(:,6) IV_MC(:,7) IV_MC(:,8) IV_MC(:,9) IV_MC(:,10)],2);

%callsDataAIR.IV_CKputs=IV_CK;
callsDataAIR.IV_MCputs2=IV_MC;
figure;
hold on;
%plot(putsDataApple.Var1, callsDataApple.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(callsDataAIR.TimeToMaturity, callsDataAIR.IV_CKputs, '-', 'DisplayName', 'IV-CK');
plot(callsDataAIR.TimeToMaturity, callsDataAIR.IVOL_MIDputs, '-', 'DisplayName', 'IVOL-MIDputs');
plot(callsDataAIR.TimeToMaturity, IV_MCmean, '-', 'DisplayName', 'IV-MCmean');
%plot(TtM, callsDataTTE.HistVol, 'DisplayName', 'Historical Volatulity of the unerlying');

title('Implied Volatility for K=149 Puts');
xlabel('Time to Maturity');
ylabel('Implied Volatility');
legend('show');
grid on;
ylim([0.15 0.35]);
hold off;

writetable(callsDataAIR, 'callsDataAIR');
callsDataAIR=readtable('callsDataAIR.txt')

%rel diff puts
relDifIVckput=abs(callsDataAIR.IV_CKputs-callsDataAIR.IVOL_MIDputs)./callsDataAIR.IVOL_MIDputs;
relDifIVmcput=abs(IV_MCmean-callsDataAIR.IVOL_MIDputs)./callsDataAIR.IVOL_MIDputs;
figure;
hold on;
plot(callsDataAIR.TimeToMaturity, relDifIVckput, '-x', 'DisplayName', 'relDifIVck');
plot(callsDataAIR.TimeToMaturity, relDifIVmcput, '-x', 'DisplayName', 'relDifIVmc');
title('Relative Difference IV');
xlabel('Time to Maturity');
ylabel('Rel Diff');
legend('show');
grid on;
ylim([0 1]);
hold off;



%calls
callsDataAIRcalls=callsDataAIR(idx:end,:);
rowsWithNaN = any(ismissing(callsDataAIRcalls), 2);
% Remove rows with NaN values
callsDataAIRcalls = callsDataAIRcalls(~rowsWithNaN, :);
callsDataAIRcalls.Calls=str2double(callsDataAIRcalls.Calls);
callsDataAIRIVcalls = readtable('Dades Opcions.xlsx','Sheet','AIRBUS').IVOL_MIDcalls;
callsDataAIRIVcalls=str2double(callsDataAIRIVcalls);
rowsWithNaN = any(ismissing(callsDataAIRIVcalls), 2);
% Remove rows with NaN values
callsDataAIRIVcalls = callsDataAIRIVcalls(~rowsWithNaN, :);

callsDataAIRcalls.IVOL_MIDcalls=callsDataAIRIVcalls;
callsDataAIRcalls.IVOL_MIDcalls=(1/100)*callsDataAIRcalls.IVOL_MIDcalls;


theoreticalPricesCKcalls=zeros(height(callsDataAIRcalls),1);
theoreticalPricesMCcalls=zeros(height(callsDataAIRcalls),1);
% Compute theoretical prices for each row
for i = 1:height(callsDataAIRcalls)
    S = callsDataAIRcalls.Underlying(i);
    T_i =callsDataAIRcalls.TimeToMaturity(i);
    K=callsDataAIRcalls.StrikeCalls(i);
    theoreticalPricesCKcalls(i) = CK_Europ(S, K, r,sigma_annualPuts,T_i, 800,1500, 'call');
    theoreticalPricesMCcalls(i) = MonteCarloEurop(S, K, r,sigma_annualPuts,T_i, 10000, 'call');
    
end
callsDataAIRcalls.theoreticalPricesCKcalls=theoreticalPricesCKcalls;
callsDataAIRcalls.theoreticalPricesMCcalls=theoreticalPricesMCcalls;


% Plot the data
figure;
hold on;
%plot(putsDataApple.Var1, callsDataApple.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(callsDataAIRcalls.Date, callsDataAIRcalls.Calls, '-x', 'DisplayName', 'Option Price');
plot(callsDataAIRcalls.Date, callsDataAIRcalls.theoreticalPricesCKcalls, '-x', 'DisplayName', 'theoreticalPricesCK');
plot(callsDataAIRcalls.Date, callsDataAIRcalls.theoreticalPricesMCcalls, '-x', 'DisplayName', 'theoreticalPricesMC');

title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
hold off;


%IV:
IV_CK=zeros(height(callsDataAIRcalls),1);
IV_MC=zeros(height(callsDataAIRcalls),10);

%IMPLIED VOL. CK
for j=1:10
for i = 1:height(callsDataAIRcalls)
    S = callsDataAIRcalls.Underlying(i);
    T_i = callsDataAIRcalls.TimeToMaturity(i);
    V=callsDataAIRcalls.Calls(i);
    K=callsDataAIRcalls.StrikeCalls(i);
    %VMC=callsDataTTE.theoreticalPricesMC(i);
    %errorFunctionCK = @(sigma) (CK_Europ(S, K, r, sigma,T_i,800,1500, 'call') - V)^2;
    errorFunctionMC = @(sigma) (MonteCarloEurop(S, K, r, sigma,T_i,100000, 'call') - V)^2;
    
    %IV_CK(i)=fminbnd(errorFunctionCK, 0.1, 1);
    IV_MC(i,j)=fminbnd(errorFunctionMC, 0.1, 1);
    
end
end
IV_MCmean=mean([IV_MC(:,1) IV_MC(:,2) IV_MC(:,3) IV_MC(:,4) IV_MC(:,5) IV_MC(:,6) IV_MC(:,7) IV_MC(:,8) IV_MC(:,9) IV_MC(:,10)],2);

%callsDataAIRcalls.IV_CKcalls=IV_CK;
callsDataAIRcalls.IV_MC=IV_MC;

figure;
hold on;
% for j =1:10
%     plot(callsDataAIRcalls.TimeToMaturity, IV_MC(:,j), '-','Color',[0.5,0.5,0.5],'HandleVisibility', 'off');
% end
plot(callsDataAIRcalls.TimeToMaturity, callsDataAIRcalls.IV_CKcalls, '-', 'DisplayName', 'IV-CK');
plot(callsDataAIRcalls.TimeToMaturity, IV_MCmean, '-', 'DisplayName', 'IV-MCmean');
plot(callsDataAIRcalls.TimeToMaturity, callsDataAIRcalls.IVOL_MIDcalls, '-', 'DisplayName', 'IVOL-MIDcalls');

title('Implied Volatility for K=164 Calls');
xlabel('Time to Maturity');
ylabel('Implied Volatility');
legend('show');
grid on;
ylim([0.14 0.28]);
hold off;

writetable(callsDataAIRcalls, 'callsDataAIRcalls');
readable('callsDataAIRcalls.txt');

relDifIVck=abs(callsDataAIRcalls.IV_CKcalls-callsDataAIRcalls.IVOL_MIDcalls)./callsDataAIRcalls.IVOL_MIDcalls;
relDifIVmc=abs(IV_MCmean-callsDataAIRcalls.IVOL_MIDcalls)./callsDataAIRcalls.IVOL_MIDcalls;
figure;
hold on;
plot(callsDataAIRcalls.TimeToMaturity, relDifIVck, '-x', 'DisplayName', 'relDifIVck');
plot(callsDataAIRcalls.TimeToMaturity, relDifIVmc, '-x', 'DisplayName', 'relDifIVmc');
title('Relative Difference IV');
xlabel('Time to Maturity');
ylabel('Rel Diff');
legend('show');
grid on;
ylim([0 1]);
hold off;