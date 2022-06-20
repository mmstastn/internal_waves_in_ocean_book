% This script makes pictures of the waves

%Pcolor plots
figure(1), clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
subplot(2,2,1)
pcolor(x,z,den),shading flat
title('Density')
subplot(2,2,2)
pcolor(x,z,vort),shading flat
title('Vorticity')
subplot(2,2,3)
pcolor(x,z,u),shading flat
title('U')
subplot(2,2,4)
pcolor(x,z,w),shading flat
title('W')

% Vertical profiles
% The filtered Ri attempts to account for regions of essentially no
% stratification.  It is ad hoc but effective.
figure(2), clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
subplot(1,2,1)
plot(uv,zv)
title('Horizontal Velocity')
subplot(1,2,2)
plot(riv,zv,'k-',rivfilt,zv,'r-')
title(['Ri (max diff ' num2str(max(abs(riv-rivfilt)),3) ')'])
legend('Raw','Filtered')
axis([0 10 min(zv) max(zv)])

