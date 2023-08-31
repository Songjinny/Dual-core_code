clc
clear
Lx=20; N=256; dx=Lx/N;
dt=0.0005; tmax=100; nmax=round(tmax/dt);
x=-Lx/2:Lx/N:Lx/2-Lx/N;dx=Lx/N;
kx=[0:N/2-1  -N/2:-1]*2*pi/Lx;
K2=kx.^2;

W0=0.5;
l=2;
% W=W0*tanh(x/l);
W=W0*sinh(x/l)./cosh(x/l).^2;


pt=0.5*x.^2+1i*W;

k=0.3;

% load('E:\dual-core\激发\data_jf\2W0=0.5_l=2_mu=0.1_p=0.5734_k=0.3_s_stable.mat','U','V');
% load('E:\dual-core\激发\data_jf\2W0=0.5_l=2_mu=-0.2_p=1.5833_k=0.3_asy_stable.mat','U','V');
% load('E:\dual-core\激发\data_jf\2W0=0.5_l=2_mu=-0.2_p=1.5349_k=0_asy_stable.mat','U','V'); %0.228---0
% load('E:\dual-core\激发\data_jf\2W0=0.5_l=2_mu=0.3_p=1.0344_k=0_s_stable.mat','U','V');
% load('E:\dual-core\激发\data_jf\2W0=0.5_l=2_mu=0.3_p=0.5232_k=0.1_asy_stable.mat','U','V');%0.2235---0
load('E:\dual-core\激发\data_jf\2W0=0.5_l=2_mu=0.1_p=0.5734_k=0.3_s_stable.mat','U','V');

u1=U;
u2=V;
% noise = max(abs(V))*(rand(1,N) + 1i*rand(1,N))*0.02;
% u1 = u1.*(1+noise);
% u2 = u2.*(1+noise);


% figure(1)
% hold on
% box on
% plot(x,real(u1),'-','linewidth',1.5);
% plot(x,imag(u1),'--','linewidth',1.5);
% plot(x,abs(u1).^2,'-.','linewidth',1.5);
% legend('Re(u)','Im(u)','|u|^2'); %axis([-7 7 -0.5 2])


udata=u1';
vdata=u2';
absuv=u1'-u2';
umax=max(abs(u1).^2); vmax=max(abs(u2).^2);
Pu=sum(abs(u1).^2)*dx; Pv=sum(abs(u2).^2)*dx;
Su=(1i/2*(u1.*ifft(1i*kx.*fft(conj(u1)))-conj(u1).*ifft(1i*kx.*fft(u1))))';
Sv=(1i/2*(u2.*ifft(1i*kx.*fft(conj(u2)))-conj(u2).*ifft(1i*kx.*fft(u2))))';
theta=abs(Pu-Pv)/(Pu+Pv);

tdata=0;

w01=0.5;
w02=0.5;
k1=0.3;
k2=0.1;

  for nn=1:nmax                               % integration begins 
   w0 = ( ((w02-w01)/2)*(1-cos(10*pi*nn*dt/tmax)) + w01).*(nn*dt<=(tmax/10))+w02.*(nn*dt>(tmax/10));  
   W=w0*sinh(x/l)./cosh(x/l).^2;
   pt=0.5*x.^2+1i*W;

%    pu=sum(abs(u1).^2)*dx; pv=sum(abs(u2).^2)*dx;
  P1=sum(abs(u1).^2)*dx; P2=sum(abs(u2).^2)*dx;
  
      k = ( ((k2-k1)/2)*(1-cos(10*pi*nn*dt/tmax)) + k1).*(nn*dt<=(tmax/10))+k2.*(nn*dt>(tmax/10)); 

    
   
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
      
    if mod(nn,round(nmax/1000)) == 0  
        udata=[udata u1']; vdata=[vdata u2'];
        absuv=[absuv u1'-u2'];
        umax=[umax, max(abs(u1).^2)];
        vmax=[vmax, max(abs(u2).^2)];
        Pu=[Pu, sum(abs(u1).^2)*dx];
        Pv=[Pv, sum(abs(u2).^2)*dx];
        Su=[Su, (1i/2*(u1.*ifft(1i*kx.*fft(conj(u1)))-conj(u1).*ifft(1i*kx.*fft(u1))))'];
        Sv=[Sv, (1i/2*(u2.*ifft(1i*kx.*fft(conj(u2)))-conj(u2).*ifft(1i*kx.*fft(u2))))'];
        theta=[theta, abs(sum(abs(u1).^2)*dx-sum(abs(u2).^2)*dx)/(sum(abs(u1).^2)*dx+sum(abs(u2).^2)*dx)];
        tdata=[tdata nn*dt];   
    end
  end  
  
  
%   
figure(1)
surf(x,tdata,(abs(udata).^2)');
colormap(jet)
shading interp
axis tight

figure(2)
surf(x,tdata,(abs(vdata).^2)');
colormap(jet)
shading interp
axis tight

% figure(3)
% plot(tdata,umax,'-','linewidth',1)
% hold on
% plot(tdata,vmax,'--','linewidth',1)

figure(4)
plot(tdata,Pu,'-','linewidth',1.5)
hold on
plot(tdata,Pv,'--','linewidth',1.5)

PP=Pu+Pv;
plot(tdata,PP,'-','linewidth',1.5)
legend('$P_1$','$P_2$','$P$'); %axis([-7 7 -0.5 2])

figure(5)
plot(tdata,theta,'-','linewidth',1.5)
% surf(x,tdata,real(Su)');
% colormap(jet)
% shading interp
% axis tight
% 
% surf(x,tdata,real(Sv)');
% colormap(jet)
% shading interp
% axis tight