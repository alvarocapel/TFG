%% ON totalEnergies
%calls
callsDataNOVO = readtable('Dades Opcions.xlsx','Sheet','NOVONORDISK','VariableNamingRule','preserve');
disp(head(callsDataNOVO));
rowsWithNaN = any(ismissing(callsDataNOVO), 2);
% Remove rows with NaN values
callsDataNOVO = callsDataNOVO(~rowsWithNaN, :);
K=900;
ExpDate = datetime('2024-06-21','InputFormat', 'yyyy-MM-dd');
r=0.02525; %German risk free rate
%volatility?
% Compute daily returns
returnsCalls = diff(log(callsDataNOVO.Underlying));
% Compute the standard deviation of the daily returns
sigma_dailyCalls = std(returnsCalls);
% Annualize the standard deviation (assuming 252 trading days in a year)
sigma_annualCalls = sigma_dailyCalls * sqrt(252);

ExpDate = zeros(height(callsDataNOVO),1)+ExpDate;

TtM=(1/252)*countBusinessDaysColumn(callsDataNOVO.Date,ExpDate);
% Initialize an array for theoretical prices
theoreticalPricesCK2 = zeros(height(callsDataNOVO), 1);
theoreticalPricesMC = zeros(height(callsDataNOVO), 1);

% Compute theoretical prices for each row
for i = 1:height(callsDataNOVO)
    S = callsDataNOVO.Underlying(i);
    T_i = TtM(i);
    theoreticalPricesCK2(i) = CK_Europ(S, K, r,sigma_annualCalls,T_i, 500,500, 'call');
    %theoreticalPricesMC(i) = MonteCarloEurop(S, K, r,sigma_annualCalls,T_i, 10000, 'call');
    
end

% Add the theoretical prices to the table
callsDataNOVO.theoreticalPricesCK = theoreticalPricesCK;
callsDataNOVO.theoreticalPricesMC = theoreticalPricesMC;
% Plot the data
figure;
hold on;
%plot(putsDataApple.Var1, callsDataApple.Underlying, '-o', 'DisplayName', 'Underlying Price');
%plot(callsDataNOVO.Date, callsDataNOVO.Calls, '-x', 'DisplayName', 'Option Price');
%plot(callsDataNOVO.Date, callsDataNOVO.theoreticalPricesCK, '-x', 'DisplayName', 'theoreticalPricesCK-r=3.5%');
%plot(callsDataNOVO.Date, callsDataNOVO.theoreticalPricesMC, '-x', 'DisplayName', 'theoreticalPricesMC');
plot(callsDataNOVO.Date,theoreticalPricesCK , '-x', 'DisplayName', 'theoreticalPricesCK-r=2.5%');
plot(callsDataNOVO.Date,theoreticalPricesCK2 , '-x', 'DisplayName', 'theoreticalPricesCK-r=2.525%');
%RHO=(theoreticalPricesCK2-theoreticalPricesCK);
%plot(callsDataNOVO.Date,RHO , '-x', 'DisplayName', 'RHO');

title('Option Market data vs Computed Black Scholes');
xlabel('Date');
ylabel('Price');
legend('show');
grid on;
%ylim([5 130]);
hold off;

callsDataNOVO.IVOL_MID=1/100*callsDataNOVO.IVOL_MID;

IV_CK=zeros(height(callsDataNOVO),1);
IV_MC=zeros(height(callsDataNOVO),10);
%IMPLIED VOL. CK
for j=1:10
tic
for i = 1:height(callsDataNOVO)
    S = callsDataNOVO.Underlying(i);
    T_i = TtM(i);
    V=callsDataNOVO.Calls(i);
    %VMC=callsDataTTE.theoreticalPricesMC(i);
    errorFunctionCK = @(sigma) (CK_Europ(S, K, r, sigma,T_i,500,500, 'call') - V)^2;
    %errorFunctionMC = @(sigma) (MonteCarloEurop(S, K, r, sigma,T_i,100000, 'call') - V)^2;
    
    IV_CK(i)=fminbnd(errorFunctionCK, 0.1, 1);
    %IV_MC(i,j)=fminbnd(errorFunctionMC, 0.05, 1.5);
    
end
timpoIVMC=toc;
end
%timpoIVCK=toc;
IV_MCmean=mean([IV_MC(:,1) IV_MC(:,2) IV_MC(:,3) IV_MC(:,4) IV_MC(:,5) IV_MC(:,6) IV_MC(:,7) IV_MC(:,8) IV_MC(:,9) IV_MC(:,10)],2);
callsDataNOVO.IV_CK=IV_CK;
callsDataNOVO.IV_MC=IV_MCmean;
callsDataNOVO.HistVol=sigma_annualCalls+zeros(height(callsDataNOVO),1);
figure;
hold on;
% for j =1:10
%     plot(TtM, IV_MC(:,j), '-','Color',[0.5,0.5,0.5],'HandleVisibility', 'off');
% end
plot(TtM, callsDataNOVO.IV_CK, '-', 'DisplayName', 'IV-CK');
plot(TtM, callsDataNOVO.IVOL_MID, '-', 'DisplayName', 'IVOL');
plot(TtM, callsDataNOVO.IV_MC, '-','DisplayName', 'IV-MCmean');

%plot(TtM, callsDataTTE.HistVol, 'DisplayName', 'Historical Volatulity of the unerlying');
title('Implied Volatility for K=900 Calls');
xlabel('Time to Maturity');
ylabel('Implied Volatility');
legend('show');
grid on;
ylim([0.1 0.6]);
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
