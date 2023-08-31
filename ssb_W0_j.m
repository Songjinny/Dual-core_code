% clc
% clear
Lx=20; N=256;
errormax=1e-12; errorCG=1e-2; c=3;
x=-Lx/2:Lx/N:Lx/2-Lx/N; dx=Lx/N;
kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;

K2=kx.^2; fftM=(c+K2);
err = 3.5*10^(-5);

l=2;
W0=[0.01:0.005:1.8];
k=[0.3];
% mu=[-0.65:0.001:0.84];
P=zeros(length(W0),length(k));
mup=zeros(length(W0),length(k));

aa=1;
i=1;
for W00=W0 
    W=W00*sinh(x/l)./cosh(x/l).^2;
%     W=W00*tanh(x/l);
    V=0.5*x.^2+1i*W;
    
    j=1;
    for kk=k
    
        U1=1.2*x.^1.*exp(-x.^2/2);
        U2=1.2*x.^1.*exp(-x.^2/2);
        mu=[aa:-0.001:3];
        for muu=mu
       
        nn=0;
        while nn<10000
           F1=V-abs(U1.*U1)-muu; G1=conj(F1);
           F2=V-abs(U2.*U2)-muu; G2=conj(F2);
           L0U1=-0.5*ifft(-K2.*fft(U1))+F1.*U1-kk*U2;
           L0U2=-0.5*ifft(-K2.*fft(U2))+F2.*U2-kk*U1;
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
             R1=-(L1A(L0U1)-kk*L0U2); R2=-(L2A(L0U2)-kk*L0U1);
             R=[R1,R2];
             MinvR1=ifft(fft(R1)./fftM); MinvR2=ifft(fft(R2)./fftM);
             MinvR=[MinvR1,MinvR2];
             R2new1=sum(conj(R1).*MinvR1); 
             R2new2=sum(conj(R2).*MinvR2);
             R201=R2new1;
             R202=R2new2;
             P1=MinvR1; P2=MinvR2;
            while((R2new1 > R201*errorCG^2) || (R2new2 > R202*errorCG^2))  
             L1P1=L1(P1)-kk*P2; L1P2=L2(P2)-kk*P1;
             LP1=L1A(L1P1)-kk*L1P2; LP2=L1A(L1P2)-kk*L1P1; 
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
            nn=nn+1;
            U1=U1+DU1;
            U2=U2+DU2;
        end
           

Nx=200;L=20;k0=2*pi/L;
G1=-2*abs(U1).^2+V-muu;
G2=-2*abs(U2).^2+V-muu;
G3=-U1.^2; G4=-U2.^2;
Gk=k*ones(size(U1));
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
            [eigvalues,I] = sort(diag(D));
             MaxImag = max(  abs( imag(eigvalues) ) );
             if MaxImag<err
                 P(i,j)=sum(abs(U1).^2)*dx+sum(abs(U2).^2)*dx; 
                 mup(i,j)=muu; 
                 aa=muu;
                 break
             end
             
            
        end
        
        j=j+1;
        
    end
    
    i=i+1;
    
    
end
% w0=W0(1:7);p=P(1:7);mmup=mup(1:7);
% save('E:\dual-core\SSB_new\data\2l=2_k=0.6.mat','W0','P')
% save('E:\dual-core\SSB_new\data\2l=2_k=0.3_s_l.mat','W0','P','mup')
% plot(k,mup(1,1:end))

plot(W0,P,'-','LineWidth',1.5)
hold on
plot(W0,p,'-','LineWidth',1.5)

plot(W0,pp,'-','LineWidth',1.5)

