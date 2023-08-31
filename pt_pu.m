clc
clear

% N = 240;            
% c = 2.2;               
% [x, D] = herdif(N, 2, c) ;   
% D2 = D(:,:,2) ; 

format long, format compact
L = 20; 
N = 512; h = 2*pi/N; 
x = h*(1:N)-h;x = L*(x-pi)/pi;
D2 = (pi/L)^2*toeplitz([-pi^2/(3*h^2)-1/6 ...
                 -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2]);
             
W0=10;
l=5;
W=W0*tanh(x/l);
% W=W0*sinh(x/l)./cosh(x/l).^2;
v=0.5*x.^2+1i*W;
k=0.1;


M=-0.5*D2+diag(v);
L_k=-k*eye(N);
MM=[M, L_k;
    L_k, M;];
NN=sparse(MM);


[V,D] = eig( M );
[M,I] = sort(diag(D));

plot(M,'.')

% u=V(:,I(2));
% plot(x,abs(u).^2)
% save('E:\dual-core\pt破缺\data\1W2.mat','u');