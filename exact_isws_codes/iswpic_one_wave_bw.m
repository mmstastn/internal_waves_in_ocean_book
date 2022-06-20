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
h1=plot(u(:,NX/2)/c,z(:,NX/2),'-');
set(h1,'Color',[0 0 0])
hold on
h2=plot(u(:,1)/c,z(:,1),'-');
set(h2,'Color',[0.8 0.8 0.8])
title('u/c')
ylabel('z (m)')
xlabel('u/c')
legend('u at crest','U(z)','Location','SouthEast')
grid on
subplot(2,2,2)
h1=plot(den(:,NX/2),z(:,NX/2),'-');
set(h1,'Color',[0 0 0])
hold on
h2=plot(den(:,1),z(:,1),'-');
set(h2,'Color',[0.8 0.8 0.8])
title('density')
ylabel('z (m)')
xlabel('\rho')
legend('\rho at crest','rhobar(z)','Location','SouthWest')
grid on
subplot(2,1,2)
h5=plot(x(end,:),u(1,:)/c,'-');
set(h5,'Color',[0 0 0])
hold on
h6=plot(x(1,:),u(end,:)/c,'-');
set(h6,'Color',[0.8 0.8 0.8])
title('u/c')
ylabel('u/c')
xlabel('x (m)')
grid on
legend('u at surface','u at bottom','Location','NorthEast')
axis([min(x(:)) max(x(:)) -0.7 0.7])
