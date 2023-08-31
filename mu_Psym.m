% clc
% clear
Lx=20;  N=256;                         % mesh parameters
max_iteration=1e5; error_tolerance=1e-12;        
x=-Lx/2:Lx/N:Lx/2-Lx/N; dx=Lx/N; kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;
k2=kx.^2; c=3;   % other parameters and initial conditions
K2=kx.^2; fftM=(c+K2);errormax=1e-12;errorCG=1e-2; 

W0=0;
l=2;
n=0;
W=W0*tanh(x/l);
% W=W0*sinh(x/l)./cosh(x/l).^2;
% if n==0
%     W=1/sqrt(2)*(6*n-3*x.*hermiteH(n,x)).*exp(-x.^2/2);
% else
%     W=1/sqrt(2)*(6*n*hermiteH(n-1,x)-3*x.*hermiteH(n,x)).*exp(-x.^2/2);
% end

V=0.5*x.^2+1i*W;


% mu=0; 

k=0.3;
% U=1.2*x.^1.*exp(-x.^2/2);                          % initial conditions
% A=sqrt( (2+sqrt(2)-4*mu)/(4*sqrt(2)) + sqrt( (2+sqrt(2)-4*mu)^2/32 - 2*k^2 ) );
% B=sqrt( (2+sqrt(2)-4*mu)/(4*sqrt(2)) - sqrt( (2+sqrt(2)-4*mu)^2/32 - 2*k^2 ) );


mu=[-2:0.1:0,0.01:0.01:0.19];
mu=[-3:0.1:-2];
P=zeros(size(mu));
theta=zeros(size(mu));
lambda=zeros(size(mu));
i=1;
for muu=mu

%     U1=1.2*exp(-x.^2/2); 
%     U2=1.2*exp(-x.^2/2);
while 1
  F1=V-abs(U1.*U1)-muu; G1=conj(F1);
  F2=V-abs(U2.*U2)-muu; G2=conj(F2);
  L0U1=-0.5*ifft(-K2.*fft(U1))+F1.*U1-k*U2;
  L0U2=-0.5*ifft(-K2.*fft(U2))+F2.*U2-k*U1;
  U=[U1,U2];
  L0U=[L0U1,L0U2];
  errorU=max(abs(L0U))/max(abs(U))
  if errorU < errormax
     break
  end
  L1= @(W) -0.5*ifft(-K2.*fft(W))+F1.*W-2*U1.*real(conj(U1).*W);
  L1A=@(W) -0.5*ifft(-K2.*fft(W))+G1.*W-2*U1.*real(conj(U1).*W);
  L2= @(W) -0.5*ifft(-K2.*fft(W))+F2.*W-2*U2.*real(conj(U2).*W);
  L2A=@(W) -0.5*ifft(-K2.*fft(W))+G2.*W-2*U2.*real(conj(U2).*W);
  
  DU1=0*U1; DU2=0*U2;
  R1=-(L1A(L0U1)-k*L0U2); R2=-(L2A(L0U2)-k*L0U1);
  R=[R1,R2];
  MinvR1=ifft(fft(R1)./fftM); MinvR2=ifft(fft(R2)./fftM);
  MinvR=[MinvR1,MinvR2];
  R2new1=sum(conj(R1).*MinvR1); 
  R2new2=sum(conj(R2).*MinvR2);
  R201=R2new1;
  R202=R2new2;
  P1=MinvR1; P2=MinvR2;
  while((R2new1 > R201*errorCG^2) || (R2new2 > R202*errorCG^2))  
    L1P1=L1(P1)-k*P2; L1P2=L2(P2)-k*P1;
    LP1=L1A(L1P1)-k*L1P2; LP2=L1A(L1P2)-k*L1P1; 
%     LP=[LP1,LP2];
    a1=R2new1/sum(real(conj(P1).*LP1));
    a2=R2new2/sum(real(conj(P2).*LP2));
    DU1=DU1+a1*P1; DU2=DU2+a2*P2;
    R1=R1-a1*LP1; R2=R2-a2*LP2; 
    MinvR1=ifft(fft(R1)./fftM); MinvR2=ifft(fft(R2)./fftM);
%     MinvR=[MinvR1,MinvR2];
    R2old1=R2new1; R2old2=R2new2;
    R2new1=sum(real(conj(R1).*MinvR1));
    R2new2=sum(real(conj(R2).*MinvR2));
    b1=R2new1/R2old1; b2=R2new2/R2old2;
    P1=MinvR1+b1*P1; P2=MinvR2+b2*P2;
  end
  
  U1=U1+DU1;
  U2=U2+DU2;
end

h = 2*pi/N;
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
D1 = (2*pi/Lx)*toeplitz(column,column([1 N:-1:2]));
column2 = [-pi^2/(3*h^2)-1/6 ...
          -0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
D2 = (2*pi/Lx)^2*toeplitz(column2); 

L1=-0.5*D2+diag(-2*abs(U1).^2+V-muu);
L2=-0.5*D2+diag(-2*abs(U2).^2+V-muu);
L3=diag(-U1.^2); L4=diag(-U2.^2);
L_k=k*eye(N); L0=0*eye(N);

M=[     L1,      -L_k,   L3,   L0;
      -L_k,        L2,   L0,   L4;
 -conj(L3),        L0,  -L1,  L_k;
        L0, -conj(L4),  L_k,  -L2;];

 [Ve,D] = eig( M);
 [eigvalues,I] = sort(diag(D));
 lambda(i) = max(  abs( imag(eigvalues) ) );
 
 P(i)=sum(abs(U1).^2)*dx+sum(abs(U2).^2)*dx;

 i=i+1;
end
hold on
box on
plot(mu(1:20),P(1:20),'-o','LineWidth',1.5)
plot(mu(20:end),P(20:end),'-b','LineWidth',1.5)
% plot(P(1:20),theta(1:20),'--b','LineWidth',1.5)
% load('E:\dual-core\SSB\data_new\0w0=0_l=2_k=0.3_mu_p_sy.mat','P','lambda','mu')

axis([-2 1 0 15])
