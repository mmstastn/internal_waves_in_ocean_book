spinsgrid2d
figure(1)
clf
betterplots
colormap temperature
for ii=1:3
subplot(3,1,ii)
ii=20+(ii-1)*20;
spinsread2dnew
contourf(x,z,u,80),shading flat,caxis([-1 1]*0.03)
hold on
contour(x,z,rho,10,'k')
end

for ii=1:3
    subplot(3,1,ii)
    x0=2.05+1.95*(ii-1);
    axis([x0-0.3 x0+0.3 0.075 0.225])
    ylabel('z (m)')
end
xlabel('x (m)')
