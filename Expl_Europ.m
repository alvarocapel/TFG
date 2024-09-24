function price= Expl_Europ(Swanted,E, r, sigma, T, N,M,type)
Smax=2*Swanted;
Smin=0;
dt=(T/(M-1)); % Time step
ds=(Smax-Smin)/(N-1); % Price step
% Black Scholes PDE with Explicit Method:

% Initializing the matrix of the option value
v=zeros(N,M);
if type == "call"
    %Initial conditions prescribed by the European Call payoff at expiry:
    % V(S,T)=max(S-E,0)
    v(1:N,1)=max((Smin+(0:N-1)*ds-E),zeros(size(1:N)))';

    % for k=1:N
    %     v(k,1)=max((Smin+(k-1)*ds-E),0)';
    % end 

    % Boundary conditions prescribed by the European Call:
    v(1,2:M)=zeros(M-1,1)'; % V(0,t)=0
    v(N,2:M)=((N-1)*ds+Smin)-E*exp(-r*((M-1)-(M-2:-1:0))*dt); % V(S,t)=S-Eexp[-r(T-t)] as S -> infininty.
elseif type == "put"
    % Initial conditions prescribed by the European Put payoff at expiry:
    % V(S, T) = max(E - S, 0)

    v(:,1) = max(E-(Smin+(0:N-1)*ds),0)';

    % Boundary conditions prescribed by the European Put:
    v(1, :) = E * exp(-r * ((M-1)-(M-1:-1:0)) * dt);   % V(S,t)=E*exp[-r(T-t)] as S -> 0.
    v(N, :) = 0; % V(S,t)=0 as S -> inf.
end


% Determining the matrix coeficients of the explicit algorithm
aa=0.5*dt*(sigma*sigma*(1:N-2).*(1:N-2)-r*(1:N-2))';
bb=1-dt*(sigma*sigma*(1:N-2).*(1:N-2)+r)';
cc=0.5*dt*(sigma*sigma*(1:N-2).*(1:N-2)+r*(1:N-2))';

% Implementing the explicit algorithm
for i=2:M
    v(2:N-1,i)=bb.*v(2:N-1,i-1)+cc.*v(3:N,i-1)+aa.*v(1:N-2,i-1); %??? (creo que ya lo pillo pero entonces me tendir que dar lo mismo.... que con los 2 fors.... no?

end

price= interp1((0:N-1)'*ds,v(:,end),Swanted);