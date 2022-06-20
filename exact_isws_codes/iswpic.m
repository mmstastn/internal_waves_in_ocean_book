figure(1+picshift)
clf
%colormap darkjet
colormap gray
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight','bold');
  umax=max(abs(u(:)));
  wmax=max(abs(w(:)));
  
subplot(2,1,1)
[conts,h]=contourf(x,z,u/umax,11);
caxis([-0.6 0.6])
set(h,'LineColor','none')
title(['u (background and wave), u max = ' num2str(umax,3) ', uamp = ' num2str(uamp)])
hold on
contour(x,z,den,6,'w')
ylabel('z (m)')
subplot(2,1,2)
[conts,h]=contourf(x,z,w/wmax,linspace(-1,1,20));
set(h,'LineColor','none')
title(['w , max w = ' num2str(wmax,3) ' uamp = ', num2str(uamp,3)])
hold on
contour(x,z,den,6,'w')
ylabel('z (m)')
xlabel('x (m)')

figure(2+picshift)
clf
%colormap darkjet
colormap gray
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight','bold');
mnrho=min(den(:));
mxrho=max(den(:));
subplot(2,2,1)
[conts,h]=contourf(x,z,den,8);
%set(h,'LineColor','none')
title('density')
ylabel('z (m)')
xlabel('x (m)')
subplot(2,2,2)
[conts,h]=contourf(x,z,ri,linspace(0,1,21));
set(h,'LineColor','none')
caxis([0.0 1])
title(['Ri, c = ' num2str(c,4)])
hold on
contour(x,z,ri,[0.249 0.251],'w')
ylabel('z (m)')
xlabel('x (m)')
denrange=max(den(:))-min(den(:));
subplot(2,2,3)
h1=plot(u(:,NX/2)/c,z(:,NX/2),'-');
set(h1,'Color',[0 0 0])
hold on
h2=plot((den(:,NX/2)-1)/denrange,z(:,NX/2),'-');
set(h2,'Color',[0.8 0.8 0.8])
ylabel('z (m)')
xlabel('scaled u and \rho')
grid on
axis([-1.1 1.1 -H 0])
legend('u','\rho')
subplot(2,2,4)
%plot(ri(:,NX/2),z(:,NX/2),'b',n2v/max(n2v),zv,'r')
h3=plot(ri(:,NX/2),z(:,NX/2),'-');
hold on
plot([0.25 0.25],[min(z(:)) max(z(:))],'k:')
h4=plot(n2v/max(n2v),zv,'Color',[0.7 0.7 0.7]);
set(h3,'Color',[0 0 0])
set(h4,'Color',[0.8 0.8 0.8])
grid on
axis([-0.05 1 -H 0])
ylabel('z (m)')
xlabel('Ri and N^2')
legend('Ri','N^2')
