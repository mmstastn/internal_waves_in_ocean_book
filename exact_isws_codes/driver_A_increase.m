%% lab scale waves using our DJL solver showing the effects of shear
clear all,close all
verbose=0; %If this is set to one you will see all the progress during solving
picshifts=[0 5 10 15];
As=[1 5 10 15]*1e-5;
figure(1)
clf
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight','bold');
colormap gray
for casei=1:length(picshifts)
    picshift=picshifts(casei); 
    uamp=0;
    A  = As(casei); % APE for wave
    NX = 512; % grid 
    NZ = 512; % grid
    H  = 0.2; %depth
    L  = 2.5; %width

    dzd=H*1e-4;dzd2=dzd*dzd;
    from_data=0;
    
    if from_data==1
     load sample_data.csv
     zd=sample_data(:,1)-20;
     dd=sample_data(:,3);
     md_density=@(z) interp1(zd,dd,z,'spline')/sample_data(1,3);
     md_d_density=@(z) (md_density(z+dzd)-md_density(z-dzd))/(2*dzd);
    else
     a_d=0.01; z0_d=0.05; d_d=0.02; % This is an example with broadening
     %a_d=0.01; z0_d=0.025; d_d=0.02; % This is an example with breaking
     %a_d=0.01; z0_d=0.025; d_d=0.005; % This is an example with low Ri
    
     md_density=@(z) 1-a_d*tanh((z+z0_d)/d_d);
     md_d_density=@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;
    end 
    zj=0.5*H;dj=0.4*H;
    md_u=@(z) uamp*tanh((z+zj)/dj);
    md_uz=@(z) (md_u(z+dzd)-md_u(z-dzd))/(2*dzd);
    md_uzz=@(z) (md_u(z+dzd)-2*md_u(z)+md_u(z-dzd))/dzd2;
    % solve problem once
    get_eta
    % %increase resolution
    x0=x; z0=z;
   
    % compute fields and save slices
    iswpost
    figure(1)
    subplot(2,2,casei)
    contourf(x,z,den,6)
    drawnow
    uths(casei,:)=u(1,:);
    ubhs(casei,:)=u(end,:);
    uvs(:,casei)=u(:,256);
    rivs(:,casei)=riv;
    cs(casei)=c;
    %iswpic_one_wave_bw
end
x1d=x(1,:);
z1d=z(:,1);
figure(1)
subplot(2,2,1)
ylabel('z (m)')
subplot(2,2,3)
ylabel('z (m)')
xlabel('x (m)')
subplot(2,2,4)
xlabel('x (m)')

figure(2)
clf
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight','bold');
for ii=1:4
    subplot(2,1,1)
    ht(ii)=plot(x1d,uths(ii,:)/cs(ii));
    hold on
    subplot(2,1,2)
    hb(ii)=plot(uvs(:,ii)/cs(ii),z1d);
    hold on
    set(ht(ii),'Color',[1 1 1]*(ii-1)*0.25)
    set(hb(ii),'Color',[1 1 1]*(ii-1)*0.25)
end
subplot(2,1,1)
axis([min(x1d) max(x1d) -0.05 1.25])
grid on
ylabel('u(x,H)/c')
xlabel('x (m)')
legend(['A = ' num2str(As(1),3)],['A = ' num2str(As(2),3)],['A = ' num2str(As(3),3)],['A = ' num2str(As(4),3)],'Location','NorthEast')
subplot(2,1,2)
axis([-1.25 1.25 min(z1d) max(z1d)])
grid on
ylabel('u(0,z)/c')
xlabel('z (m)')
legend(['A = ' num2str(As(1),3)],['A = ' num2str(As(2),3)],['A = ' num2str(As(3),3)],['A = ' num2str(As(4),3)],'Location','SouthEast')

figure(3)
clf
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight','bold');
for ii=1:4
    subplot(2,1,1)
    hb(ii)=plot(uvs(:,ii)/cs(ii),z1d);
    hold on
    set(hb(ii),'Color',[1 1 1]*(ii-1)*0.25)
    subplot(2,1,2)
    h3(ii)=plot(rivs(:,ii),z1d);
    hold on
    set(h3(ii),'Color',[1 1 1]*(ii-1)*0.25)
end
subplot(2,1,1)
axis([-1.25 1.25 min(z1d) max(z1d)])
grid on
ylabel('u(L/2,z)/c')
xlabel('z (m)')
legend(['A = ' num2str(As(1),3)],['A = ' num2str(As(2),3)],['A = ' num2str(As(3),3)],['A = ' num2str(As(4),3)],'Location','SouthEast')
subplot(2,1,2)
axis([0 1 min(z1d) max(z1d)])
grid on
ylabel('Ri(L/2,z)/c')
xlabel('z (m)')
legend(['A = ' num2str(As(1),3)],['A = ' num2str(As(2),3)],['A = ' num2str(As(3),3)],['A = ' num2str(As(4),3)],'Location','SouthEast')
