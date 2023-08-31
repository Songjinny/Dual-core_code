clc
clear
Lx=40; N=256;
dt=0.0005; tmax=1000; nmax=round(tmax/dt);
x=-Lx/2:Lx/N:Lx/2-Lx/N;dx=Lx/N;
kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;
K2=kx.^2;

W0=0.1;
l=2;

W=W0*sinh(x/l)./cosh(x/l).^2;


pt=0*x.^2+1i*W;

k=0.1;


% load('E:\dual-core\V=0\W0=0.1_l=2_mu=-1_k=0.mat','U','V');
% load('E:\dual-core\V=0\data\aW0=0.1_l=2_mu=-1_k=0.1_asy_unstable.mat','U','V');

% load('E:\dual-core\V=0\data\W0=0.1_l=2_mu=-0.5_k=0.1_asy_stable.mat','U','V'); %稳定
% load('E:\dual-core\V=0\data\W0=0.1_l=2_mu=-0.3_k=0.1_asy_stable.mat','U','V');
% load('E:\dual-core\V=0\data\aW0=0.1_l=2_mu=-0.5_k=0.1_asy_stable.mat','U','V');
% load('E:\dual-core\V=0\data\aW0=0.2_l=2_mu=-0.5_k=0.1_asy_stable.mat','U','V');

load('E:\dual-core\V=0\data\aW0=0.1_l=2_mu=-0.5_k=0.1_asy_stable.mat','U','V');
load('E:\dual-core\V=0\data\W0=0.1_l=2_mu=-0.18_k=0.1_sy_stable.mat','U','V');

u1=U;
u2=V;
noise = (rand(1,N) + 1i*rand(1,N))*0;
u1 = u1.*(1+noise);
u2 = u2.*(1+noise);

% u1 = u1+0.01*(f1+conj(g1));
% u2 = u2+0.01*(f2+conj(g2));

% figure(1)
% hold on
% box on
% plot(x,real(u1),'-','linewidth',1.5);
% plot(x,imag(u1),'--','linewidth',1.5);
% plot(x,abs(u1).^2,'-.','linewidth',1.5);
% legend('Re(u)','Im(u)','|u|^2'); %axis([-7 7 -0.5 2])


udata=u1';  vdata=u2';
umax=max(abs(u1).^2); vmax=max(abs(u2).^2);
tdata=0;


  for nn=1:nmax                               % integration begins 
    du1=-1i*(-0.5*ifft(-K2.*fft(u1))+pt.*u1-abs(u1).^2.*u1-k*u2);  
    dv1=-1i*(-0.5*ifft(-K2.*fft(u2))+pt.*u2-abs(u2).^2.*u2-k*u1); 
    v1=u1+0.5*du1*dt; v2=u2+0.5*dv1*dt; 
    du2=-1i*(-0.5*ifft(-K2.*fft(v1))+pt.*v1-abs(v1).^2.*v1-k*v2);  
    dv2=-1i*(-0.5*ifft(-K2.*fft(v2))+pt.*v2-abs(v2).^2.*v2-k*v1); 
    v1=u1+0.5*du2*dt; v2=u2+0.5*dv2*dt; 
    du3=-1i*(-0.5*ifft(-K2.*fft(v1))+pt.*v1-abs(v1).^2.*v1-k*v2);  
    dv3=-1i*(-0.5*ifft(-K2.*fft(v2))+pt.*v2-abs(v2).^2.*v2-k*v1); 
    v1=u1+    du3*dt; v2=u2+    dv3*dt; 
    du4=-1i*(-0.5*ifft(-K2.*fft(v1))+pt.*v1-abs(v1).^2.*v1-k*v2); 
    dv4=-1i*(-0.5*ifft(-K2.*fft(v2))+pt.*v2-abs(v2).^2.*v2-k*v1); 
    u1=u1+(du1+2*du2+2*du3+du4)*dt/6;
    u2=u2+(dv1+2*dv2+2*dv3+dv4)*dt/6;
      
    if mod(nn,round(nmax/300)) == 0  
        udata=[udata u1']; vdata=[vdata u2']; 
        umax=[umax, max(abs(u1).^2)];
        vmax=[vmax, max(abs(u2).^2)];
        tdata=[tdata nn*dt];   
    end
  end   
figure(1)
surf(x,tdata(1:151),(abs(udata(:,1:151)).^2)');
colormap(jet)
shading interp
axis tight

figure(2)
surf(x,tdata(1:151),(abs(vdata(:,1:151)).^2)');
colormap(jet)
shading interp
axis tight

plot(tdata,umax,'-','linewidth',1)
hold on
plot(tdata,vmax,'-','linewidth',1)