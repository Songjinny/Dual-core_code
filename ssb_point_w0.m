% clc
% clear
Lx=20; N=256;
errormax=1e-12; errorCG=1e-2; c=3;
x=-Lx/2:Lx/N:Lx/2-Lx/N; dx=Lx/N;
kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;

K2=kx.^2; fftM=(c+K2);
err = 3.5*10^(-8);


W0=2.5;
l=[2];

k=[0.1:0.01:1];
P=zeros(length(l),length(k));
mup=zeros(length(l),length(k));

aa=0.95;
i=1;
for ll=l 
    W=W0*sinh(x/ll)./cosh(x/ll).^2;
%     W=W0*tanh(x/ll);
    V=0.5*x.^2+1i*W;
    
    j=1;
    for kk=k
        
        mu=[aa:-0.001:-2];
        
        for muu=mu
        
%         U1=1.2*exp(-x.^2/2);
%         U2=1.2*exp(-x.^2/2);
        
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
             L1= @(S) -0.5*ifft(-K2.*fft(S))+F1.*S-2*U1.*real(conj(U1).*S);
             L1A=@(S) -0.5*ifft(-K2.*fft(S))+G1.*S-2*U1.*real(conj(U1).*S);
             L2= @(S) -0.5*ifft(-K2.*fft(S))+F2.*S-2*U2.*real(conj(U2).*S);
             L2A=@(S) -0.5*ifft(-K2.*fft(S))+G2.*S-2*U2.*real(conj(U2).*S);
  
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
           

            h = 2*pi/N;
            column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
            D1 = (2*pi/Lx)*toeplitz(column,column([1 N:-1:2]));
            column2 = [-pi^2/(3*h^2)-1/6 ...
                         -0.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
            D2 = (2*pi/Lx)^2*toeplitz(column2); 

            L1=-0.5*D2+diag(-2*abs(U1).^2+V-muu);
            L2=-0.5*D2+diag(-2*abs(U2).^2+V-muu);
            L3=diag(-U1.^2); L4=diag(-U2.^2);
            L_k=kk*eye(N); L0=0*eye(N);

            M=[     L1,      -L_k,   L3,   L0;
                 -L_k,        L2,   L0,   L4;
            -conj(L3),        L0,  -L1,  L_k;
                   L0, -conj(L4),  L_k,  -L2;];

            [Ve,D] = eig( M);
            [eigvalues,I] = sort(diag(D));
             MaxImag = max(  abs( imag(eigvalues) ) );
             if MaxImag>err
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

plot(k,P,'--','LineWidth',1.5)
hold on
plot(k,mup,'--','LineWidth',1.5)
% save('E:\dual-core\SSB\data_new\2w0=2.5_l=2.mat','P','k','mup')

% view(90,0);
% W=W0*tanh(x/l);

% if n==0
%     W=1/sqrt(2)*(6*n-3*x.*hermiteH(n,x)).*exp(-x.^2/2);
% else
%     W=1/sqrt(2)*(6*n*hermiteH(n-1,x)-3*x.*hermiteH(n,x)).*exp(-x.^2/2);
% end


