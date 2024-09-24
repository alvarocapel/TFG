%% ON totalEnergies
%calls
callsDataTTE = readtable('Dades Opcions.xlsx','Sheet','TOTALENERGIES');
disp(callsDataTTE);
rowsWithNaN = any(ismissing(callsDataTTE), 2);
% Remove rows with NaN values
callsDataTTE = callsDataTTE(~rowsWithNaN, :);
K=68;
ExpDate = datetime('2024-06-21','InputFormat', 'yyyy-MM-dd');
r=0.025; %German risk free rate
%volatility?
% Compute daily returns
returnsCalls = diff(log(callsDataTTE.Underlying));
% Compute the standard deviation of the daily returns
sigma_dailyCalls = std(returnsCalls);
% Annualize the standard deviation (assuming 252 trading days in a year)
sigma_annualCalls = sigma_dailyCalls * sqrt(252);

ExpDate = zeros(height(callsDataTTE),1)+ExpDate;

TtM=(1/252)*countBusinessDaysColumn(callsDataTTE.Date,ExpDate);
% Initialize an array for theoretical prices
theoreticalPricesCK = zeros(height(callsDataTTE), 1);
theoreticalPricesMC = zeros(height(callsDataTTE), 1);

% Compute theoretical prices for each row
for i = 1:height(callsDataTTE)
    S = callsDataTTE.Underlying(i);
    T_i = TtM(i);
    theoreticalPricesCK(i) = CK_Europ(S, K, r,sigma_annualCalls,T_i, 800,1500, 'call');
    theoreticalPricesMC(i) = MonteCarloEurop(S, K, r,sigma_annualCalls,T_i, 10000, 'call');
    
end

% Add the theoretical prices to the table
callsDataTTE.theoreticalPricesCK = theoreticalPricesCK;
callsDataTTE.theoreticalPricesMC = theoreticalPricesMC;
% Plot the data
figure;
hold on;
plot(callsDataTTE.Date,zeros(height(callsDataTTE),1)+68 , '-', 'DisplayName', 'Strike');
plot(callsDataTTE.Date, callsDataTTE.Underlying, '-x', 'DisplayName', 'Underlying Price');
% plot(callsDataTTE.Date, callsDataTTE.theoreticalPricesCK, '-x', 'DisplayName', 'theoreticalPricesCK');
% plot(callsDataTTE.Date, callsDataTTE.theoreticalPricesMC, '-x', 'DisplayName', 'theoreticalPricesMC');

title('Underlying Price');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
hold off;

callsDataTTE.IVOL_MIDcalls=1/100*callsDataTTE.IVOL_MIDcalls;
callsDataTTE.IVOL_MIDputs=1/100*callsDataTTE.IVOL_MIDputs;

IV_CK=zeros(height(callsDataTTE),1);
IV_MC=zeros(height(callsDataTTE),10);
%IMPLIED VOL. CK
tic
for j=1:10

for i = 1:height(callsDataTTE)
    S = callsDataTTE.Underlying(i);
    T_i = TtM(i);
    V=callsDataTTE.Call(i);
    %VMC=callsDataTTE.theoreticalPricesMC(i);
    %errorFunctionCK = @(sigma) (CK_Europ(S, K, r, sigma,T_i,800,1500, 'call') - V)^2;
    errorFunctionMC = @(sigma) (MonteCarloEurop(S, K, r, sigma,T_i,100000, 'call') - V)^2;
    
    %IV_CK(i)=fminbnd(errorFunctionCK, 0.05, 1.5);
    IV_MC(i,j)=fminbnd(errorFunctionMC, 0.05, 1.5);
    
end

end
timpoIVMC=toc;
%timpoIVCK=toc;
IV_MCmean=mean([IV_MC(:,1) IV_MC(:,2) IV_MC(:,3) IV_MC(:,4) IV_MC(:,5) IV_MC(:,6) IV_MC(:,7) IV_MC(:,8) IV_MC(:,9) IV_MC(:,10)],2);
%callsDataTTE.IV_CK=IV_CK;
callsDataTTE.IV_MC=IV_MC;
callsDataTTE.HistVol=sigma_annualCalls+zeros(height(callsDataTTE),1);
figure;
hold on;
% for j =1:10
%     plot(TtM, IV_MC(:,j), '-','Color',[0.5,0.5,0.5],'HandleVisibility', 'off');
% end
plot(TtM, callsDataTTE.IV_CK, '-', 'DisplayName', 'IV-CK');
plot(TtM, callsDataTTE.IVOL_MIDcalls, '-', 'DisplayName', 'IVOL');
plot(TtM, IV_MCmean, '-', 'Color', 'y','DisplayName', 'IV-MCmean');

%plot(TtM, callsDataTTE.HistVol, 'DisplayName', 'Historical Volatulity of the unerlying');
title('Implied Volatility for K=68 Calls');
xlabel('Time to Maturity');
ylabel('Implied Volatility');
legend('show');
grid on;
ylim([0.12 0.26]);
hold off;

callsDataTTE.IV_MCmean=IV_MCmean;
relDifIVck=abs(callsDataTTE.IV_CK-callsDataTTE.IVOL_MIDcalls)./callsDataTTE.IVOL_MIDcalls;
relDifIVmc=abs(callsDataTTE.IV_MCmean-callsDataTTE.IVOL_MIDcalls)./callsDataTTE.IVOL_MIDcalls;
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

%puts
K=65;
ExpDate = datetime('2024-06-21','InputFormat', 'yyyy-MM-dd');
r=0.025; %German risk free rate
%volatility?
% Compute daily returns
returnsCalls = diff(log(callsDataTTE.Underlying));
% Compute the standard deviation of the daily returns
sigma_dailyCalls = std(returnsCalls);
% Annualize the standard deviation (assuming 252 trading days in a year)
sigma_annualCalls = sigma_dailyCalls * sqrt(252);

ExpDate = zeros(height(callsDataTTE),1)+ExpDate;

TtM=(1/252)*countBusinessDaysColumn(callsDataTTE.Date,ExpDate);

theoreticalPricesCK = zeros(height(callsDataTTE), 1);
theoreticalPricesMC = zeros(height(callsDataTTE), 1);

% Compute theoretical prices for each row
for i = 1:height(callsDataTTE)
    S = callsDataTTE.Underlying(i);
    T_i = TtM(i);
    theoreticalPricesCK(i) = CK_Europ(S, K, r,sigma_annualCalls,T_i, 600,600, 'put');
    theoreticalPricesMC(i) = MonteCarloEurop(S, K, r,sigma_annualCalls,T_i, 10000, 'put');
    
end

% Add the theoretical prices to the table
callsDataTTE.theoreticalPricesCKputs = theoreticalPricesCK;
callsDataTTE.theoreticalPricesMCputs = theoreticalPricesMC;

% Plot the data
figure;
hold on;
%plot(putsDataApple.Var1, callsDataApple.Underlying, '-o', 'DisplayName', 'Underlying Price');
plot(callsDataTTE.Date, callsDataTTE.Put, '-x', 'DisplayName', 'Option Price');
plot(callsDataTTE.Date, callsDataTTE.theoreticalPricesCKputs, '-x', 'DisplayName', 'theoreticalPricesCK');
plot(callsDataTTE.Date, callsDataTTE.theoreticalPricesMCputs, '-x', 'DisplayName', 'theoreticalPricesMC');

title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
hold off;
plot(callsDataTTE.Date, MSE, 'x', 'DisplayName', 'errors ^2');
plot(callsDataTTE.Date, mean(MSE)+zeros(size(callsDataTTE.Date)), '-', 'DisplayName', 'MSE');


for i = 1:height(callsDataTTE)
    
    CK=callsDataTTE.theoreticalPricesCKputs(i);
    real=callsDataTTE.Put(i);
    MSE(i)=(CK-real)^2;
end
mean(MSE);

%REAL IV?
% callsDataTTEIVcalls = readtable('Dades Opcions.xlsx','Sheet','TOTALENERGIES').IVOL_MIDcalls;
% rowsWithNaN = any(ismissing(callsDataTTEIVcalls), 2);
% callsDataTTEIVputs = readtable('Dades Opcions.xlsx','Sheet','TOTALENERGIES').IVOL_MIDputs;
% rowsWithNaN = any(ismissing(callsDataTTEIVputs), 2);
% % Remove rows with NaN values
% callsDataTTEIVcalls = callsDataTTEIVcalls(~rowsWithNaN, :);
% callsDataTTEIVputs = callsDataTTEIVputs(~rowsWithNaN, :);
% 
% callsDataTTE.IVOL_MIDcalls=callsDataTTEIVcalls;
% callsDataTTE.IVOL_MIDputs=callsDataTTEIVputs;
% callsDataTTE.IVOL_MIDcalls=(1/100)*callsDataTTE.IVOL_MIDcalls;
% callsDataTTE.IVOL_MIDputs=(1/100)*callsDataTTE.IVOL_MIDputs;

IV_CKputs=zeros(height(callsDataTTE),1);
IV_MC=zeros(height(callsDataTTE),10);
%IMPLIED VOL. CK
for j=1:10
tic
for i = 1:height(callsDataTTE)
    S = callsDataTTE.Underlying(i);
    T_i = TtM(i);
    V=callsDataTTE.Put(i);
    %VMC=callsDataTTE.theoreticalPricesMC(i);
    %errorFunctionCK = @(sigma) (CK_Europ(S, K, r, sigma,T_i,800,800, 'put') - V)^2;
    errorFunctionMC = @(sigma) (MonteCarloEurop(S, K, r, sigma,T_i,10000, 'put') - V)^2;
    
    %IV_CK(i)=fminbnd(errorFunctionCK, 0.05, 1.5);
    IV_MC(i,j)=fminbnd(errorFunctionMC, 0.05, 1.5);
    
end
tiempoIVCKputs=toc;
end
callsDataTTE.IV_CKputs=IV_CK;

IV_MCmean=mean([IV_MC(:,1) IV_MC(:,2) IV_MC(:,3) IV_MC(:,4) IV_MC(:,5) IV_MC(:,6) IV_MC(:,7) IV_MC(:,8) IV_MC(:,9) IV_MC(:,10)],2);

callsDataTTE.IV_MCputs=IV_MCmean;
%rel diff puts
relDifIVckput=abs(callsDataTTE.IV_CKputs-callsDataTTE.IVOL_MIDputs)./callsDataTTE.IVOL_MIDputs;
relDifIVmcput=abs(callsDataTTE.IV_MCputs-callsDataTTE.IVOL_MIDputs)./callsDataTTE.IVOL_MIDputs;
figure;
hold on;
plot(TtM, callsDataTTE.IVOL_MIDputs, '-x', 'DisplayName', 'relDifIVck');
plot(TtM, callsDataTTE.IVOL_MIDcalls, '-x', 'DisplayName', 'relDifIVmc');
title('Relative Difference IV');
xlabel('Time to Maturity');
ylabel('Rel Diff');
legend('show');
grid on;
%ylim([0 1]);
hold off;
