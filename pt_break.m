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
             
err = 10^(-9);

lmin = 0.1;   lmax = 2;   lin = lmax - lmin;
W0min = 0 ;   W0max = 10;   W0in = W0max - W0min;
   
k=0.5;
Ml = 100;
MW0 = 100;
dl = lin / Ml;
dW0 = W0in / MW0;

l = (lmin : dl : lmax);
W0 = (W0min : dW0 : W0max);

Ml = Ml+1;
MW0 = MW0+1;

map0s = zeros(1,Ml);

top = 1;
parfor jj = 1 : Ml 

   top = 1 ;    
   signn=1; 
    while signn==1&&top<=MW0

%         W=W0(top)*tanh(x/l(jj));
      W=W0(top)*sinh(x/l(jj))./cosh(x/l(jj)).^2;

        u = 0.5*x.^2+1i*W;
        M=-0.5*D2+diag(u);
        
%         L_k=-k*eye(N);
%         MM=[M, L_k;
%             L_k, M;];
        % calculate the eigenvalues and judge it
        [V,D] = eig( M );
        if max(  abs( imag(diag(D)) ) ) < err
            top = top + 1;
        else
            signn = 0;
        end
    end
    
    
    map0s(jj) = top-1;
end
% for each oemgaPT(jj), the boundary is dw0*top
hold on;
plot(l,dW0*map0s','-','linewidth',1);
xlabel('l'); ylabel('W_{0}');
plot(l,-dW0*map0s','-','linewidth',1);
axis([lmin 2 0  W0max])
box on
% save('E:\dual-core\pt破缺\data\W2.mat','dW0','map0s','l')
% y1 = dW0*map0s;
% y2 = - dW0*map0s;
% fill([V1,fliplr(V1)],[y1,fliplr(y2)],[0.667,0.75,1])             