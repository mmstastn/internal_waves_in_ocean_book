%% This code plots KdV solitary waves and relevant parameters for two layer
%% theory as the layer depeths and wave amplitudes are changed.

g=9.81;
delrho=0.01; %percent density change
H=100;
h1s=linspace(0.51*H,0.99*H,91); % bottom layer thickness
B0s=linspace(0.01*H,0.25*H,91);% amplitde, 
[h1h1 B0B0]=meshgrid(h1s,B0s);
h2h2=H-h1h1;
c0c0=sqrt(g*delrho*h1h1.*h2h2/H);
alpha=1.5*c0c0.*(h1h1-h2h2)./(h1h1.*h2h2);
beta=(1/6)*c0c0*h1h1*h2h2;
cwnl=c0c0+(B0B0/3).*alpha;
gamma=sqrt(B0B0.*alpha./(12*beta));
widthwnl=1./gamma;

figure(1)
clf
subplot(2,1,1)
h=plot(B0s/H,cwnl(:,10)./c0c0(:,10),'k-',B0s/H,cwnl(:,55)./c0c0(:,55),'k--',B0s/H,cwnl(:,90)./c0c0(:,90),'k:')
set(h,'linewidth',2)
xlabel('B_0/H')
ylabel('c/c_0')
h1a=h1s(10)/H;
h1b=h1s(55)/H;
h1c=h1s(90)/H;
grid on
legend(['h_1/H = ' num2str(h1a,3)],['h_1/H = ' num2str(h1b,3)],['h_1/H = ' num2str(h1c,3)],'Location','Northwest')
subplot(2,1,2)
h=plot(h1s/H,cwnl(10,:)./c0c0(10,:),'k-',h1s/H,cwnl(55,:)./c0c0(55,:),'k--',h1s/H,cwnl(:,90)./c0c0(:,90),'k:')
set(h,'linewidth',2)
xlabel('h_1/H')
ylabel('c/c_0')
B0a=B0s(10)/H;
B0b=B0s(55)/H;
B0c=B0s(90)/H;
grid on
legend(['B_0/H = ' num2str(B0a,3)],['B_0/H = ' num2str(B0b,3)],['B_0/H = ' num2str(B0c,3)],'Location','Northwest')
axis([0.5 1 0 4])

figure(2)
clf
subplot(2,1,1)
h=plot(B0s/H,widthwnl(:,10)./H,'k-',B0s/H,widthwnl(:,55)./H,'k--',B0s/H,widthwnl(:,90)./H,'k:')
set(h,'linewidth',2)
xlabel('B_0/H')
ylabel('wave width / H')
h1a=h1s(10)/H;
h1b=h1s(55)/H;
h1c=h1s(90)/H;
grid on
legend(['h_1/H = ' num2str(h1a,3)],['h_1/H = ' num2str(h1b,3)],['h_1/H = ' num2str(h1c,3)],'Location','Northeast')
subplot(2,1,2)
h=plot(h1s/H,widthwnl(10,:)./H,'k-',h1s/H,widthwnl(55,:)./H,'k--',h1s/H,widthwnl(:,90)./H,'k:')
set(h,'linewidth',2)
xlabel('h_1/H')
ylabel('wave width / H')
B0a=B0s(10)/H;
B0b=B0s(55)/H;
B0c=B0s(90)/H;
grid on
legend(['B_0/H = ' num2str(B0a,3)],['B_0/H = ' num2str(B0b,3)],['B_0/H = ' num2str(B0c,3)],'Location','Northeast')
%axis([0.5 1 0 4])

myx=linspace(-1e4,1e4,1001);
figure(3)
clf
subplot(2,1,1)
h=plot(myx/H,(B0s(90)/H)*sech(gamma(90,85)*myx).^2,'k-',myx/H,(B0s(55)/H)*sech(gamma(55,85)*myx).^2,'k--',myx/H,(B0s(10)/H)*sech(gamma(10,85)*myx).^2,'k:')
set(h,'linewidth',2)
xlabel('x/H')
ylabel('solitary wave profile')
title(['h_1/H = ' num2str(h1s(85)/H,4)])
legend(['B_0/H = ' num2str(B0s(90)/H,4)],['B_0/H = ' num2str(B0s(55)/H,4)],['B_0/H = ' num2str(B0s(10)/H,4)])
grid on
subplot(2,1,2)
h=plot(myx/H,(B0s(55)/H)*sech(gamma(55,85)*myx).^2,'k-',myx/H,(B0s(55)/H)*sech(gamma(55,65)*myx).^2,'k--')
set(h,'linewidth',2)
xlabel('x/H')
ylabel('solitary wave profile')
title(['B_0/H = ' num2str(B0s(55)/H,4)])
legend(['h_1/H = ' num2str(h1s(85)/H,4)],['h_1/H = ' num2str(h1s(65)/H,4)])
grid on