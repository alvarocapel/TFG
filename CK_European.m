function price= CK_European(Swanted,E, r, sigma, T, M,N,type)
Smax=2*Swanted;
Smin=0;
dt=(T/(M-1)); % Time step
ds=(Smax-Smin)/(N-1); % Price step

% Initializing the matrix of the option value
v(1:N,1:M) = 0.0;
v=zeros(N,M);
if type == "call"
    %Initial conditions of European Call payoff at expiry:V(S,T)=max(S-E,0)
    v(1:N,1)=max((Smin+(0:N-1)*ds-E),0)';

    % Boundary conditions prescribed by the European Call:
    v(1,2:M)=zeros(M-1,1)'; % V(0,t)=0
    v(N,2:M)=((N-1)*ds+Smin)-E*exp(-r*((M-1)-(M-2:-1:0))*dt); 
    %V(S,t)=S-Eexp[-r(T-t)] as S -> inf
    
elseif type == "put"
    %Initial conditions prescribed by the European Put payoff at expiry:
    % V(S,T)=max(E-S,0)
    v(:,1)=max(E-(Smin+(0:N-1)*ds),zeros(size(1:N)))';

    % Boundary conditions prescribed by the European PUT:
    v(1, :) = E * exp(-r * (((M-1)-(M-1:-1:0)) * dt)); %V(S,t)=E*exp[-r(T-t)] as S -> 0.
    v(N,:)= 0;
end

for n=1:N-1 %price. From n=0 to n= 1 before last node
    alphaa(n)=0.5*0.5*dt*(r*(n-1)-(n-1)^2*sigma^2);
    betaa(n)=1+0.5*dt*((n-1)^2*sigma^2+r);
    gammaa(n)=-0.5*0.5*dt*(r*(n-1)+(n-1)^2*sigma^2);
end
% termes independents:
for n=1:N-1 %price
    b1(n)=0.5*0.5*(sigma^2*(n-1)^2-r*(n-1))*dt;
    b2(n)=1-0.5*(sigma^2*(n-1)^2+r)*dt;
    b3(n)=0.5*0.5*(sigma^2*(n-1)^2+r*(n-1))*dt;
end

% Constructing the tridiagonal matrix
A1 = diag(betaa(2:end)) + diag(alphaa(3:end), -1) +diag(gammaa(2:end-1),1);
[L,U]=lu(A1);
aux=zeros(N-2,1);
aux2=zeros(N-2,1);
B = diag(b2(2:end)) + diag(b1(3:end), -1) + diag(b3(2:end-1), 1);
for j = 2:M  %now the backward behavior is implemented in the for loop
    
        aux(1)= b1(2)*v(1,j); %1st unknown node coef * boundary node
        aux(end)=b3(end)*v(end,j); %last unknown node coef * boundary node  
        aux2(1)=b1(2)*v(1,j-1);
        aux2(end)=b3(end)*v(end,j-1);
        
        % Solving the system of equations
        v(2:end-1, j) = U\(L\(B*v(2:end-1,j-1))+aux+aux2);
        
end
price= interp1(Smin+(0:N-1)'*ds,v(:,end),Swanted);
end 