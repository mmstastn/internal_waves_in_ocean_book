%%FFT based method for the BBM equation
%% B_t + c_0 B_x + alpha B B_x - (beta/c)B_xxt = 0
% Uses leapfrog scheme

clear all; close all
%Physical parameters
H=100;
h1=20; %upper layer depth
h2=H-h1; %lower layer depth
g=9.81; drho=0.005; gp=g*drho;
ctwolayer=sqrt(gp*h1*h2/H);
betatwolayer=ctwolayer*h1*h2/6;
alphatwolayer=1.5*ctwolayer*(h1-h2)/(h1*h2);

%define a spatial grid
xmin = -1e4;
xmax = 1e4;
N = 2^9;
x = linspace(xmin,xmax,N+1); x=x(1:end-1);
dx = x(2)-x(1);

%make initial condition
B1 = -0.1*H*sech(x/(0.025*xmax)).^2;
B2 = -B1;%.*sin(2*pi*x/2);
u0=B1;


%make wave numbers
nyquist_freq = 2*pi/(xmax-xmin);
ks=[0:N/2-1 0 -N/2+1:-1]*nyquist_freq;
ks2=ks.*ks; ks3=ks2.*ks;

% define the filter; this section is presently unused
% Here I use hypervisocsity so the order should be even
% The parameter can be somewhat automated but here I tune by hand
% filtord=8;
% hypervisc=1e-16;
% myfilt=1./(1+hypervisc*ks.^filtord);

% figure(100)
% clf
% subplot(2,1,1)
% plot(ks,'bo')
% ylabel('wavenumber')
% subplot(2,1,2)
% plot(fftshift(ks),fftshift(myfilt),'bo')
% xlabel('wavenumber')
% ylabel('filter')

%time step and number of steps
t=0;
dt = 1e-1;twodt=2*dt;
numstps=1000;
numouts=200;
figure(1)
clf
% this sets the thick line plotting parameters useful for hard copies and journals
 set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
subplot(2,1,1)
plot(x,u0,'k-')
subplot(2,1,2)
u0f=fft(u0);spec0=log10(u0f.*conj(u0f));
plot(fftshift(ks),fftshift(spec0),'k.-')

bbmfact=1./(1+(betatwolayer/ctwolayer).*ks2);
%start with one backward Euler time step then move to leapfrog
B1lin=-ctwolayer*sqrt(-1)*ks.*fft(B1);
B1nl=-0.5*alphatwolayer*sqrt(-1)*ks.*fft(B1.^2);
B1p=B1-dt*real(ifft(bbmfact.*(B1lin+B1nl)));
B2lin=-ctwolayer*sqrt(-1)*ks.*fft(B2);
B2nl=-0.5*alphatwolayer*sqrt(-1)*ks.*fft(B2.^2);
B2p=B2-dt*real(ifft(bbmfact.*(B2lin+B2nl)));
for ii=1:numouts
   for jj=1:numstps 
    t=t+dt;
    % explicit method for BBM
    B1lin=-ctwolayer*sqrt(-1)*ks.*fft(B1);
    B1nl=-0.5*alphatwolayer*sqrt(-1)*ks.*fft(B1.^2);
    B1f = B1p+twodt*real(ifft(bbmfact.*(B1lin+B1nl)));
    % explicit method for BBM
    B2lin=-ctwolayer*sqrt(-1)*ks.*fft(B2);
    B2nl=-0.5*alphatwolayer*sqrt(-1)*ks.*fft(B2.^2);
    B2f = B2p+twodt*real(ifft(bbmfact.*(B2lin+B2nl)));
    % now rotate the times
    B1p=B1; B1=B1f; B2p=B2; B2=B2f;
   end 
   figure(2)
   clf
% this sets the thick line plotting parameters useful for hard copies and journals
% old Matlab only
%  set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
%         'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
%           'DefaultAxesFontWeight','bold');
   %subplot(2,1,1)
   h=plot(x/H,u0/H,'k:',x/H,B1/H,'k-',x/H,B2/H,'k--')
   set(h,'linewidth',2)
   grid on
   xlabel('x/H');
   ylabel('B /H');
   title(['time = ' num2str(t,2)]);
   axis([xmin/H xmax/H -1.5*max(max(abs(u0)))/H 1.5*max(max(abs(u0)))/H])
   drawnow

end
