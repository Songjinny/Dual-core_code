% clc
% clear
Lx=20;  N=256;                         % mesh parameters
max_iteration=1e5; error_tolerance=1e-12;        
x=-Lx/2:Lx/N:Lx/2-Lx/N; dx=Lx/N; kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;
k2=kx.^2;

W0=0;
l=2;
n=0;
% W=W0*tanh(x/l);
W=W0*sinh(x/l)./cosh(x/l).^2;
% if n==0
%     W=1/sqrt(2)*(6*n-3*x.*hermiteH(n,x)).*exp(-x.^2/2);
% else
%     W=1/sqrt(2)*(6*n*hermiteH(n-1,x)-3*x.*hermiteH(n,x)).*exp(-x.^2/2);
% end

pt=0.5*x.^2+1i*W;


% mu=0; 
c=3.8; DT=0.01;   % other parameters and initial conditions
k=0.2;
% U=1.2*x.^1.*exp(-x.^2/2);                          % initial conditions
% A=sqrt( (2+sqrt(2)-4*mu)/(4*sqrt(2)) + sqrt( (2+sqrt(2)-4*mu)^2/32 - 2*k^2 ) );
% B=sqrt( (2+sqrt(2)-4*mu)/(4*sqrt(2)) - sqrt( (2+sqrt(2)-4*mu)^2/32 - 2*k^2 ) );
% U=1.2*x.^2.*exp(-x.^2/2).*exp(1i*sqrt(pi/2)*erf(x/2)); 
% U=2*x.^4.*exp(-x.^2/1).*exp(1i*sqrt(pi/2)*erf(x/2));
% V=1.2*x.^2.*exp(-x.^2/2).*exp(1i*sqrt(pi/2)*erf(x/2));

mu=[-7:0.1:0,0.02:0.02:0.08];
P=zeros(size(mu));
theta=zeros(size(mu));
lambda=zeros(size(mu));
i=1;
for muu=mu

for nn=1:max_iteration                             % MSOM iterations start                   
    Uold=U;
    Vold=V;
    L0U=-0.5*ifft(-k2.*fft(U))+(pt-abs(U).^2-muu).*U-k*V; 
    L0V=-0.5*ifft(-k2.*fft(V))+(pt-abs(V).^2-muu).*V-k*U; 
    MinvL0U=ifft(fft(L0U)./(k2+c)); 
    MinvL0V=ifft(fft(L0V)./(k2+c)); 
    L1HermitMinvL0U=-0.5*ifft(-k2.*fft(MinvL0U))+(conj(pt)-abs(U).^2-muu).*MinvL0U-2*U.*real(conj(U).*MinvL0U)-k*MinvL0V; 
    L1HermitMinvL0V=-0.5*ifft(-k2.*fft(MinvL0V))+(conj(pt)-abs(V).^2-muu).*MinvL0V-2*V.*real(conj(V).*MinvL0V)-k*MinvL0U; 
    MinvL1HermitMinvL0U=ifft(fft(L1HermitMinvL0U)./(k2+c));
    MinvL1HermitMinvL0V=ifft(fft(L1HermitMinvL0V)./(k2+c));

    if nn == 1
        U=U-MinvL1HermitMinvL0U*DT; 
        V=V-MinvL1HermitMinvL0V*DT; 
    else
        L1G1=-0.5*ifft(-k2.*fft(G1))+(pt-abs(U).^2-muu).*G1-2*U.*real(conj(U).*G1)-k*G2;
        L1G2=-0.5*ifft(-k2.*fft(G2))+(pt-abs(V).^2-muu).*G2-2*V.*real(conj(V).*G2)-k*G1;
        MinvL1G1=ifft(fft(L1G1)./(k2+c)); 
        MinvL1G2=ifft(fft(L1G2)./(k2+c)); 
        MG1=ifft(fft(G1).*(k2+c)); 
        MG2=ifft(fft(G2).*(k2+c)); 
        alpha1=1/sum(conj(MG1).*G1)-1/(DT*sum(conj(L1G1).*MinvL1G1)); 
        alpha2=1/sum(conj(MG2).*G2)-1/(DT*sum(conj(L1G2).*MinvL1G2));
        innerproduct1=sum(real(G1).*real(L1HermitMinvL0U)+imag(G1).*imag(L1HermitMinvL0U)); 
        innerproduct2=sum(real(G2).*real(L1HermitMinvL0V)+imag(G2).*imag(L1HermitMinvL0V)); 
        U=U-(MinvL1HermitMinvL0U-alpha1*innerproduct1*G1)*DT; 
        V=V-(MinvL1HermitMinvL0V-alpha2*innerproduct2*G2)*DT; 
    end
    
    G1=U-Uold; 
    G2=V-Vold;
    Uerror(nn)=sqrt(sum(abs(U-Uold).^2+abs(V-Vold).^2)*dx); Uerror(nn)
    if Uerror(nn) < error_tolerance 
         break
    end
end

h = 2*pi/N;
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
D1 = (2*pi/Lx)*toeplitz(column,column([1 N:-1:2]));
column2 = [-pi^2/(3*h^2)-1/6 ...
          -0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
D2 = (2*pi/Lx)^2*toeplitz(column2); 

L1=-0.5*D2+diag(-2*abs(U).^2+pt-muu);
L2=-0.5*D2+diag(-2*abs(V).^2+pt-muu);
L3=diag(-U.^2); L4=diag(-V.^2);
L_k=k*eye(N); L0=0*eye(N);

M=[     L1,      -L_k,   L3,   L0;
      -L_k,        L2,   L0,   L4;
 -conj(L3),        L0,  -L1,  L_k;
        L0, -conj(L4),  L_k,  -L2;];

 [Ve,D] = eig( M);
 [eigvalues,I] = sort(diag(D));
 lambda(i) = max(  abs( imag(eigvalues) ) );
 
 P(i)=sum(abs(U).^2)*dx+sum(abs(V).^2)*dx;
 theta(i)=abs(sum(abs(U).^2)*dx-sum(abs(V).^2)*dx)/(sum(abs(U).^2)*dx+sum(abs(V).^2)*dx);
 i=i+1;
end



plot(P,theta,'-b','LineWidth',1.5)
% load('E:\dual-core\SSB\data_theta\2w0=1_l=2_k=0.3.mat','P','lambda','theta')
hold on
box on
plot(P(1:11),theta(1:11),'--r','LineWidth',1.5)
plot(P(11:end),theta(11:end),'-b','LineWidth',1.5)

