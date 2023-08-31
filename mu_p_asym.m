% clc
% clear
Lx=40;  N=512;                         % mesh parameters
max_iteration=1e5; error_tolerance=1e-12;        
x=-Lx/2:Lx/N:Lx/2-Lx/N; dx=Lx/N; kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;
k2=kx.^2;

W0=0.1;
l=2;

% W=W0*tanh(x/l);
W=W0*sinh(x/l)./cosh(x/l).^2;
% if n==0
%     W=1/sqrt(2)*(6*n-3*x.*hermiteH(n,x)).*exp(-x.^2/2);
% else
%     W=1/sqrt(2)*(6*n*hermiteH(n-1,x)-3*x.*hermiteH(n,x)).*exp(-x.^2/2);
% end

pt=0*x.^2+1i*W;


% mu=0; 
c=3.8; DT=0.01;   % other parameters and initial conditions
k=0.1;
% U=1.2*x.^1.*exp(-x.^2/2);                          % initial conditions
% A=sqrt( (2+sqrt(2)-4*mu)/(4*sqrt(2)) + sqrt( (2+sqrt(2)-4*mu)^2/32 - 2*k^2 ) );
% B=sqrt( (2+sqrt(2)-4*mu)/(4*sqrt(2)) - sqrt( (2+sqrt(2)-4*mu)^2/32 - 2*k^2 ) );
% U=1.2*exp(-x.^2/2).*exp(1i*sqrt(pi/2)*erf(x/2)); 
% U=1.2*x.^1.*exp(-x.^2/1).*exp(-1i*sqrt(pi/2)*erf(x/2));
% V=1*x.^0.*exp(-x.^2/2).*exp(-1i*sqrt(pi/2)*erf(x/2));

mu=[-2:0.1:-0.3,-0.29:0.01:-0.16];
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

Nx=200;L=40;k0=2*pi/L;
G1=-2*abs(U).^2+pt-muu;
G2=-2*abs(V).^2+pt-muu;
G3=-U.^2; G4=-V.^2;
Gk=k*ones(size(U));
for n=-Nx:1:Nx
c1(n+Nx+1)=dx*sum(G1.*exp(-1i*k0*n*x))/L;
c2(n+Nx+1)=dx*sum(G2.*exp(-1i*k0*n*x))/L;
c3(n+Nx+1)=dx*sum(G3.*exp(-1i*k0*n*x))/L;
c4(n+Nx+1)=dx*sum(G4.*exp(-1i*k0*n*x))/L;
ck(n+Nx+1)=dx*sum(Gk.*exp(-1i*k0*n*x))/L;
end
D=(1i*k0*diag([-Nx:Nx])).^2;
C0=zeros(2*Nx+1);
C1=toeplitz([c1(Nx+1:2*Nx+1) zeros(1,Nx)],[c1(Nx+1:-1:1) zeros(1,Nx)]);
C2=toeplitz([c2(Nx+1:2*Nx+1) zeros(1,Nx)],[c2(Nx+1:-1:1) zeros(1,Nx)]);
C3=toeplitz([c3(Nx+1:2*Nx+1) zeros(1,Nx)],[c3(Nx+1:-1:1) zeros(1,Nx)]);
C4=toeplitz([c4(Nx+1:2*Nx+1) zeros(1,Nx)],[c4(Nx+1:-1:1) zeros(1,Nx)]);
Ck=toeplitz([ck(Nx+1:2*Nx+1) zeros(1,Nx)],[ck(Nx+1:-1:1) zeros(1,Nx)]);
M=[     -0.5*D+C1,      -Ck,   C3,   C0;
      -Ck,        -0.5*D+C2,   C0,   C4;
 -conj(C3),        C0,  -(-0.5*D+C1),  Ck;
        C0, -conj(C4),  Ck,  -(-0.5*D+C2);];
 [Ve,D] = eig( M);
%  [Ve,D] = eigs( M,500,'SM');
 [eigvalues,I] = sort(diag(D));
 lambda(i) = max(  abs( imag(eigvalues) ) );
 
 P(i)=sum(abs(U).^2)*dx+sum(abs(V).^2)*dx;
 theta(i)=abs(sum(abs(U).^2)*dx-sum(abs(V).^2)*dx)/(sum(abs(U).^2)*dx+sum(abs(V).^2)*dx);
 i=i+1;
end

% load('E:\dual-core\V=0\data_mup\w0=0.1_l=2_k=0.1_mu_p_asy.mat','P','lambda','mu')
% load('E:\dual-core\V=0\data_mup\w0=0.1_l=2_k=0.1_mu_p_asy1.mat','P','theta','lambda','mu')








