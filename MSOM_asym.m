% clc
% clear
Lx=20;  N=256;                         % mesh parameters
max_iteration=5e5; error_tolerance=1e-12;        
x=-Lx/2:Lx/N:Lx/2-Lx/N; dx=Lx/N; kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;
k2=kx.^2;

W0=0;
l=2;
% W=W0*tanh(x/l);
W=W0*sinh(x/l)./cosh(x/l).^2;


pt=0.5*x.^2+1i*W;

c=3.8; DT=0.01;   % other parameters and initial conditions
mu=-10; 

k=0.3;
% U=1.2*x.^1.*exp(-x.^2/2);                          % initial conditions

% % % U=W0*x.^0.*exp(-x.^2/2).*exp(1i*sqrt(pi/2)*erf(x/2)); 
% U=1.2*x.^1.*exp(-x.^2/1).*exp(-1i*sqrt(pi/2)*erf(x/2));
% V=1*x.^0.*exp(-x.^2/2).*exp(-1i*sqrt(pi/2)*erf(x/2));


for nn=1:max_iteration                             % MSOM iterations start                   
    Uold=U;
    Vold=V;
    L0U=-0.5*ifft(-k2.*fft(U))+(pt-abs(U).^2-mu).*U-k*V; 
    L0V=-0.5*ifft(-k2.*fft(V))+(pt-abs(V).^2-mu).*V-k*U; 
    MinvL0U=ifft(fft(L0U)./(k2+c)); 
    MinvL0V=ifft(fft(L0V)./(k2+c)); 
    L1HermitMinvL0U=-0.5*ifft(-k2.*fft(MinvL0U))+(conj(pt)-abs(U).^2-mu).*MinvL0U-2*U.*real(conj(U).*MinvL0U)-k*MinvL0V; 
    L1HermitMinvL0V=-0.5*ifft(-k2.*fft(MinvL0V))+(conj(pt)-abs(V).^2-mu).*MinvL0V-2*V.*real(conj(V).*MinvL0V)-k*MinvL0U; 
    MinvL1HermitMinvL0U=ifft(fft(L1HermitMinvL0U)./(k2+c));
    MinvL1HermitMinvL0V=ifft(fft(L1HermitMinvL0V)./(k2+c));

    if nn == 1
        U=U-MinvL1HermitMinvL0U*DT; 
        V=V-MinvL1HermitMinvL0V*DT; 
    else
        L1G1=-0.5*ifft(-k2.*fft(G1))+(pt-abs(U).^2-mu).*G1-2*U.*real(conj(U).*G1)-k*G2;
        L1G2=-0.5*ifft(-k2.*fft(G2))+(pt-abs(V).^2-mu).*G2-2*V.*real(conj(V).*G2)-k*G1;
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

L1=-0.5*D2+diag(-2*abs(U).^2+pt-mu);
L2=-0.5*D2+diag(-2*abs(V).^2+pt-mu);
L3=diag(-U.^2); L4=diag(-V.^2);
L_k=k*eye(N); L0=0*eye(N);

M=[     L1,      -L_k,   L3,   L0;
      -L_k,        L2,   L0,   L4;
 -conj(L3),        L0,  -L1,  L_k;
        L0, -conj(L4),  L_k,  -L2;];

 [Ve,D] = eig( M);
%  [Ve,D] = eigs( M,100,'SM');
 [eigvalues,I] = sort(diag(D));
 MaxImag = max(  abs( imag(eigvalues) ) )
% figure(1) 
% plot(eigvalues,'o') 
 
 
 %--------Fourier collacation method
 
% Nx=140;L=20;k0=2*pi/L;
% G1=-2*abs(U).^2+pt-mu;
% G2=-2*abs(V).^2+pt-mu;
% G3=-U.^2; G4=-V.^2;
% Gk=k*ones(size(U));
% for n=-Nx:1:Nx
% c1(n+Nx+1)=dx*sum(G1.*exp(-1i*k0*n*x))/L;
% c2(n+Nx+1)=dx*sum(G2.*exp(-1i*k0*n*x))/L;
% c3(n+Nx+1)=dx*sum(G3.*exp(-1i*k0*n*x))/L;
% c4(n+Nx+1)=dx*sum(G4.*exp(-1i*k0*n*x))/L;
% ck(n+Nx+1)=dx*sum(Gk.*exp(-1i*k0*n*x))/L;
% end
% D=(1i*k0*diag([-Nx:Nx])).^2;
% C0=zeros(2*Nx+1);
% C1=toeplitz([c1(Nx+1:2*Nx+1) zeros(1,Nx)],[c1(Nx+1:-1:1) zeros(1,Nx)]);
% C2=toeplitz([c2(Nx+1:2*Nx+1) zeros(1,Nx)],[c2(Nx+1:-1:1) zeros(1,Nx)]);
% C3=toeplitz([c3(Nx+1:2*Nx+1) zeros(1,Nx)],[c3(Nx+1:-1:1) zeros(1,Nx)]);
% C4=toeplitz([c4(Nx+1:2*Nx+1) zeros(1,Nx)],[c4(Nx+1:-1:1) zeros(1,Nx)]);
% Ck=toeplitz([ck(Nx+1:2*Nx+1) zeros(1,Nx)],[ck(Nx+1:-1:1) zeros(1,Nx)]);
% M=[     -0.5*D+C1,      -Ck,   C3,   C0;
%       -Ck,        -0.5*D+C2,   C0,   C4;
%  -conj(C3),        C0,  -(-0.5*D+C1),  Ck;
%         C0, -conj(C4),  Ck,  -(-0.5*D+C2);];
%  [Ve,D] = eig( M);
% %  [Ve,D] = eigs( M,500,'SM');
%  [eigvalues,I] = sort(diag(D));
%  MaxImag = max(  abs( imag(eigvalues) ) )
%  figure(2) 
% plot(eigvalues,'o') 
 
 
 
 




figure(1)                                           % plotting results
hold on
box on
plot(x,abs(U).^2,'-','linewidth',1.5)
plot(x,abs(V).^2,'--','linewidth',1.5)
% plot(x,real(U),'-','linewidth',1.5)
% plot(x,real(V),'-','linewidth',1.5)
% plot(x,imag(U),'-.','linewidth',1.5)
% plot(x,imag(V),'-.','linewidth',1.5)
% legend('Re$(u)$','Re$(v)$','Im$(u)$','Im$(v)$'); %axis([-7 7 -0.5 2])
legend('$|u|^2$','$|v|^2$'); %axis([-7 7 -0.5 2])
axes('position',[0.15 0.6 0.3 0.3]) 
box on
hold on
plot(eigvalues,'o') 
% axis([-10 10 -1 1])


per=Ve(:,3)';
f1=per(1:N);
f2=per(N+1:2*N);
g1=per(2*N+1:3*N);
g2=per(3*N+1:end);
e1=eigvalues(3,1);

% save('E:\dual-core\演化\data_j\2W0=0.2_l=2_mu=0.7_p=2.4818_k=0.3_asy_unstable.mat','U','V');


% save('E:\dual-core\演化\data\1W0=2.8_l=2_mu=0.4_p=1.32_k=0.3_asy_unstable_p.mat','f1','f2','g1','g2','e1');

% figure(3)            
% semilogy(1:length(Uerror), Uerror, 'linewidth', 2)

L0U=-0.5*ifft(-k2.*fft(U))+(pt-abs(U).^2-mu).*U-k*V; 
L0V=-0.5*ifft(-k2.*fft(V))+(pt-abs(V).^2-mu).*V-k*U; 
error=max(abs(L0U))+max(abs(L0V))

P=sum(abs(U).^2)*dx+sum(abs(V).^2)*dx
theta=abs((sum(abs(U).^2)*dx-sum(abs(V).^2)*dx))/P


% save('E:\dual-core\演化\data\2W0=2.5_l=1_mu=-0.5_k=0.5.mat','U','V');


% W0=0.5;
% l=1;
% k=0.15;
% mu=0;
% chi=sqrt(1+W0/l);
% L=2*W0^2*sqrt(chi)/(pi*l*sqrt(l^2*chi+2))+(chi^2+1)/(4*chi);
% AA=sqrt(  ((L-mu)+sqrt((L-mu)^2-4*k^2)) /sqrt(2)  );
% BB=sqrt(  ((L-mu)-sqrt((L-mu)^2-4*k^2)) /sqrt(2)  );
% uu=AA*exp(-chi/2*x.^2).*exp(-1i*W0*erf(x/l));
% vv=BB*exp(-chi/2*x.^2).*exp(-1i*W0*erf(x/l));
% plot(x,abs(uu).^2,'--','linewidth',1.5,'color',[0.85,0.33,0.10])
% hold on
% plot(x,abs(vv).^2,'--','linewidth',1.5,'color',[0.85,0.33,0.10])
% 
% axes('position',[0.15 0.6 0.3 0.3]) 
% box on
% 
% plot(x,real(U),'linewidth',1.5)
% hold on
% plot(x,imag(U),'linewidth',1.5,'color',[0.93,0.69,0.13])
% hold on
% plot(x,real(vv),'--','linewidth',1.5,'color',[0.85,0.33,0.10])
% plot(x,imag(vv),'--','linewidth',1.5,'color',[0.49,0.18,0.56])
% plot(x,abs(uu).^2,'--','linewidth',1.5,'color',[0.85,0.33,0.10])