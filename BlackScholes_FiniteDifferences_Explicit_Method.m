%% Black Scholes PDE with Explicit Method:

% Parameters of the problem:
%clear all
r=0.2; % Interest rate
sigma=0.25; % Volatility of the underlying
M=1600+1; % Number of time points
N=160+1; % Number of share price points
Smax=20; % Maximum share price considered
Smin=0; % Minimum share price considered
T=1; % Maturation (expiry)of contract
E=10; % Exercise price of the underlying

dt=(T/(M-1)); % Time step
ds=(Smax-Smin)/(N-1); % Price step

% Initializing the matrix of the option value
v(1:N,1:M) = 0.0;
v=zeros(N,M);

%Initial conditions prescribed by the European Call payoff at expiry:
% V(S,T)=max(S-E,0)
v(1:N,1)=max((Smin+(0:N-1)*ds-E),zeros(size(1:N)))';

% for k=1:N
%     v(k,1)=max((Smin+(k-1)*ds-E),0)';
% end 

% Boundary conditions prescribed by the European Call:
v(1,2:M)=zeros(M-1,1)'; % V(0,t)=0
v(N,2:M)=((N-1)*ds+Smin)-E*exp(-r*((M-1)-(M-2:-1:0))*dt); % V(S,t)=S-Eexp[-r(T-t)] as S -> infininty.

%v1 code
v1=v;
%time
for n=1:N %price
    aaa(n)=0.5*(sigma^2*(n-1)^2-r*(n-1))*dt;
    bbb(n)=1-(sigma^2*(n-1)^2+r)*dt;
    ccc(n)=0.5*(sigma^2*(n-1)^2+r*(n-1))*dt;
end


for j=2:M
    for n=2:N-1
        
        v1(n,j)=aaa(n)*v1(n-1,j-1)+bbb(n)*v1(n,j-1)+ccc(n)*v1(n+1,j-1);  %explicit method,
        
    end
end
v1=fliplr(v1);
Swanted=Smax*0.5;
price= interp1((0:N-1)'*ds,v1(:,1),Swanted)
%% v code

% Determining the matrix coeficients of the explicit algorithm
aa=0.5*dt*(sigma*sigma*(1:N-2).*(1:N-2)-r*(1:N-2))';
bb=1-dt*(sigma*sigma*(1:N-2).*(1:N-2)+r)';
cc=0.5*dt*(sigma*sigma*(1:N-2).*(1:N-2)+r*(1:N-2))';

% Implementing the explicit algorithm
for i=2:M
    v(2:N-1,i)=bb.*v(2:N-1,i-1)+cc.*v(3:N,i-1)+aa.*v(1:N-2,i-1); %??? (creo que ya lo pillo pero entonces me tendir que dar lo mismo.... que con los 2 fors.... no?

end
spy(v-v1)


% Reversal of the time components in the matrix as the solution of the BlackScholes
% equation was performed backwards
v=fliplr(v);
% Figure of the value of the option, V(S,t), as a function of S
% at three different times:t=0, T/2 and T (expiry).

figure(1)
plot(Smin+ds*(0:N-1),v(1:N,1)','r-',Smin+ds*(0:N-1),v(1:N,round(M/2))','g-',Smin+ds*(0:N-1),v(1:N,M)','b-');
xlabel('S');
ylabel('V(S,t)');
title('European Call Option within the Explicit Method');
% Figure of the Value of the option, V(S,t)
figure(2)
mesh(Smin+ds*(0:N-1),dt*(0:M-1),v(1:N,1:M)')
title('European Call Option value, V(S,t), within the Explicit Method')
xlabel('S')
ylabel('t')

%% My code for Explicit Put Option:
% Parameters of the problem:
r = 0.2; % Interest rate
sigma = 0.25; % Volatility of the underlying
M = 1600+1; % Number of time points
N = 160+1; % Number of share price points
Smax = 20; % Maximum share price considered
Smin = 0; % Minimum share price considered
T = 1; % Maturation (expiry) of contract
E = 10; % Exercise price of the underlying

dt = T /(M-1); % Time step
ds = (Smax - Smin) /(N-1); % Price step

% Initializing the matrix of the option value
v = zeros(N, M);

% Initial conditions prescribed by the European Put payoff at expiry:
% V(S, T) = max(E - S, 0)

v(:,1) = max(E-(Smin+(0:N-1)*ds),0)';

% Boundary conditions prescribed by the European Put:
v(1, :) = E * exp(-r * ((M-1:-1:0)) * dt);   % V(S,t)=E*exp[-r(T-t)] as S -> 0.
v(N, :) = 0; % V(S,t)=0 as S -> inf.

%v1 code
v1 = v;

for n=1:N %price
    aaa(n)=0.5*(sigma^2*(n-1)^2-r*(n-1))*dt;
    bbb(n)=1-(sigma^2*(n-1)^2+r)*dt;
    ccc(n)=0.5*(sigma^2*(n-1)^2+r*(n-1))*dt;

end
% veti=0:N;
% a = 0.5*dt*(sigma^2*veti - r).*veti;
% b = 1- dt*(sigma^2*veti.^2 + r);
% c = 0.5*dt*(sigma^2*veti + r).*veti;
for j = 2:M
    for n = 2:N-1
        v1(n, j) =  aaa(n) * v1(n - 1, j - 1) + bbb(n) * v1(n, j - 1) + ccc(n) * v1(n + 1, j - 1);  %explicit method
    end
end
v1=fliplr(v1);
Swanted=Smax*0.5;
price= interp1((0:N-1)'*ds,v1(:,1),Swanted)

% Determining the matrix coeficients of the explicit algorithm
aa=0.5*dt*(sigma*sigma*(1:N-1).*(1:N-1)-r*(1:N-1))';
bb=1-dt*(sigma*sigma*(1:N-1).*(1:N-1)+r)';
cc=0.5*dt*(sigma*sigma*(1:N-1).*(1:N-1)+r*(1:N-1))';

% Implementing the explicit algorithm
for i=2:M+1
    v(2:N,i)=bb.*v(2:N,i-1)+cc.*v(3:N+1,i-1)+aa.*v(1:N-1,i-1); %??? (creo que ya lo pillo pero entonces me tendir que dar lo mismo.... que con los 2 fors.... no?

end


