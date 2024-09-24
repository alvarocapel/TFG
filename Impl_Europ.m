function price= Impl_Europ(Swanted,E, r, sigma, T, N,M,type)
Smax=2*Swanted;
Smin=0;
dt=(T/(M-1)); % Time step
ds=(Smax-Smin)/(N-1); % Price step


% Initializing the matrix of the option value
v(1:N,1:M) = 0.0;
v=zeros(N,M);
if type == "call" 
    %Initial conditions prescribed by the European Call payoff at expiry:
    % V(S,T)=max(S-E,0)
    v(:,1)=max((Smin+(0:N-1)*ds-E),zeros(size(1:N)))';

    % for k=1:N
    %     v(k,1)=max((Smin+(k-1)*ds-E),0)';
    % end 

    v(1,:)=0; % V(0,t)=0
    v(N,:)=((N-1)*ds+Smin)-E*exp(-r*((M-1)-(M-1:-1:0))*dt); % V(S,t)=S-Eexp[-r(T-t)] as S -> infinty.

elseif type == "put"
    %Initial conditions prescribed by the European Put payoff at expiry:
    % V(S,T)=max(E-S,0)
    v(:,1)=max(E-(Smin+(0:N-1)*ds),zeros(size(1:N)))';

    % for k=1:N
    %     v(k,1)=max((Smin+(k-1)*ds-E),0)';
    % end 


    % Boundary conditions prescribed by the European PUT:
    v(1, :) = E * exp(-r * (((M-1)-(M-1:-1:0)) * dt)); % V(S,t)=E*exp[-r(T-t)] as S -> 0.
    v(N,:)= 0;
end

for n=1:N-1 %price. From n=1 to n=last node
    alphaa(n)=-0.5*dt*(sigma^2*(n-1)^2-r*(n-1));
    betaa(n)=1-(-(dt*(sigma^2*(n-1)^2+r)));
    gammaa(n)=-0.5*dt*(sigma^2*(n-1)^2+r*(n-1));
end
% Implementing the implicit algorithm

% Constructing the tridiagonal matrix
A1 = diag(betaa(2:end)) + diag(alphaa(3:end), -1) + diag(gammaa(2:end-1), 1);
[L,U]=lu(A1);
aux=zeros(N-2,1);

for j = 2:M  %now the backward behavior is implemented in the for loop
    aux(1)= -alphaa(2)*v(1,j); % el coeficiente del 1st unknown node * el valor del boundary node
    aux(end)=-gammaa(end)*v(end,j); %coef del last unknown node * valor boundary node

    % Solving the system of equations
    v(2:end-1, j) = U\(L\ (v(2:end-1, j-1)+aux));
    
end
%v1=fliplr(v1);
Swanted=Smax*0.5;
price= interp1((0:N-1)'*ds,v(:,end),Swanted);

end