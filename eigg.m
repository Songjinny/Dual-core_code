clc
clear


format long, format compact
L = 20; 
N = 512; h = 2*pi/N; 
x = h*(1:N)-h;x = L*(x-pi)/pi;
D2 = (pi/L)^2*toeplitz([-pi^2/(3*h^2)-1/6 ...
                 -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2]);  




l=2;
% V0min = -15;   V0max = 2;   V0in = V0max - V0min;
W0min = 0 ;   W0max = 4.5;   W0in = W0max - W0min;
  

% MV0 = 10  ;
MW0 = 150 ;
% dV0 = V0in / MV0;
dW0 = W0in / MW0;

% searching meshes
% V0 = (V0min : dV0 : V0max);
W0 = (W0min : dW0 : W0max);

% N intervals --> N+1 points
% MV0 = MV0+1;
MW0 = MW0+1;

% map0s: record the broken/unbroken boundary
Eigs = [];


top = 1;
for jj = 1 : MW0 

       
%         w=W0(jj)*tanh(x/l);
        w=W0(jj)*sinh(x/l)./cosh(x/l).^2;

        u = 0.5*x.^2 + 1i*w;
        % calculate the eigenvalues and judge it
        Eigs = [Eigs,sort(eig( -0.5*D2+diag(u) ))] ;
        
end
    

%     real
figure(1)
hold on
    plot(W0,real(Eigs(1,:)),'LineWidth',2);
    xlabel('W_{0}'); ylabel('Re(\lambda)');
    plot(W0,real(Eigs(2,:)),'--','LineWidth',2);
    plot(W0,real(Eigs(3,:)),'-.','LineWidth',2);
    plot(W0,real(Eigs(4,:)),'.','LineWidth',2);
box on
% figure(2)
  %imag
    eig = imag(Eigs(1:2,:));
    temp=0;
    for jj=1:MW0-1
        if eig(1,jj)>eig(2,jj)
           temp= eig(1,jj);
           eig(1,jj)=eig(2,jj);
           eig(2,jj)=temp;
        end
    end
    
    eig2 = imag(Eigs(3:4,:));
    temp=0;
    for jj=1:MW0-1
        if eig2(1,jj)>eig2(2,jj)
           temp= eig2(1,jj);
           eig2(1,jj)=eig2(2,jj);
           eig2(2,jj)=temp;
        end
    end
    
    hold on;box on;
   
    plot(W0,eig(1,:),'LineWidth',2);
    xlabel('W_{0}'); ylabel('Im(\lambda)');
    plot(W0,eig(2,:),'--','LineWidth',2);
      plot(W0,eig2(1,:),'-.','LineWidth',2);
        plot(W0,eig2(2,:),'.','LineWidth',2);



