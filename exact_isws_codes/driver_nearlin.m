%% lab scale waves using our DJL solver showing the behaviour
%% near the linear stratification for which the DJL linearizes
clear all,close all
verbose=0;

NX = 256; % grid 
NZ = 256; % grid
H  = 0.2; %depth
L  = 2; %width

start_time = clock;
dzd=H*1e-4;dzd2=dzd*dzd;
uamp=0.0;
from_data=0;
dels=0.009:-0.001:0.001;
umids=zeros(256,9);
usurfs=zeros(9,256);

targamp=0.15;eps=0.001;
for deli=1:9  
    A  = 1e-5; % APE for wave
if from_data==1
 load sample_data.csv
 zd=sample_data(:,1)-20;
 dd=sample_data(:,3);
 md_density=@(z) interp1(zd,dd,z,'spline')/sample_data(1,3);
 md_d_density=@(z) (md_density(z+dzd)-md_density(z-dzd))/(2*dzd);
else
 delrho=0.02; a_d=dels(deli); z0_d=0.05; d_d=0.01;  n00=delrho/H;c00=sqrt(n00)*H/pi;
% zpyc=-0.05;dpyc=0.01;drho_halved=0.02;
 md_density=@(z) 1-delrho*z/H-a_d*tanh((z+z0_d)/d_d);
 md_d_density=@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;
end 
zj=0.5*H;dj=0.4*H;
md_u=@(z) uamp*tanh((z+zj)/dj);
md_uz=@(z) (md_u(z+dzd)-md_u(z-dzd))/(2*dzd);
md_uzz=@(z) (md_u(z+dzd)-2*md_u(z)+md_u(z-dzd))/dzd2;
% solve problem once for low and once for high
get_eta
x0=x;z0=z;
ampnow=max(abs(eta(:)))/H
while abs(ampnow-targamp)>eps
    A=A/(ampnow/targamp);
    get_eta
    ampnow=max(abs(eta(:)))/H
end
iswpost
etamax(deli)=max(abs(eta(:)));
cs(deli)=c;
umids(:,deli)=u(:,128);
usurfs(deli,:)=u(1,:);
end

end_time=clock;
fprintf('Total wall clock time: %f\n',etime(end_time, start_time));
figure(1)
clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight','bold');
subplot(2,2,1)
plot(umids/c00,z(:,128)/H)
hold on
plot(umids(:,end)/c00,z(:,128)/H,'k--','linewidth',2)
ylabel('z')
xlabel('scaled u(0,z)')
axis([-2 2 -1 0])
subplot(2,2,2)
plot([dels/0.02 0],[cs c00]/c00,'bo-')
ylabel('scaled prop. speed')
xlabel('scaled pycnocline strength')
axis([0 0.5 1 5])
subplot(2,1,2)
plot(x(1,:)/H,usurfs/c00)
hold on
plot(x(1,:)/H,usurfs(end,:)/c00,'k--','linewidth',2)
axis([0 L/H 0 2])
xlabel('x')
ylabel('scaled u(x,0)')


