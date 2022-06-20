%% lab scale waves using our DJL solver showing the effects of shear
clear all,close all
verbose=0; %If this is set to one you will see all the progress during solving
uamps=[-0.05 0 0.05];
figure(1)
clf
colormap gray
for casei=1:length(uamps)
    uamp=uamps(casei);
    A  = 4e-5; % APE for wave
    NX = 512; % grid 
    NZ = 512; % grid
    H  = 0.2; %depth
    L  = 3; %width

    dzd=H*1e-4;dzd2=dzd*dzd;
    from_data=0;
    if from_data==1
     load sample_data.csv
     zd=sample_data(:,1)-20;
     dd=sample_data(:,3);
     md_density=@(z) interp1(zd,dd,z,'spline')/sample_data(1,3);
     md_d_density=@(z) (md_density(z+dzd)-md_density(z-dzd))/(2*dzd);
    else
     a_d=0.01; z0_d=0.05; d_d=0.02;  
     md_density=@(z) 1-a_d*tanh((z+z0_d)/d_d);
     md_d_density=@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;
    end 
    zj=0.5*H;dj=0.4*H;
    md_u=@(z) uamp*tanh((z+zj)/dj);
    md_uz=@(z) (md_u(z+dzd)-md_u(z-dzd))/(2*dzd);
    md_uzz=@(z) (md_u(z+dzd)-2*md_u(z)+md_u(z-dzd))/dzd2;
    % solve problem once
    get_eta
    iswpost
    x0=x;z0=z;
    % plot things
    figure(1)
    subplot(2,3,casei+3)
    contourf(x,z,u/c,10)
    caxis([-1 1]*0.6)
    subplot(2,3,casei)
    contourf(x,z,den,6)
end
subplot(2,3,1)
ylabel('z (m)')
title('Negative shear')
subplot(2,3,2)
title('No shear')
subplot(2,3,3)
title('Positive shear')
subplot(2,3,4)
ylabel('z (m)')
xlabel('x (m)')
subplot(2,3,5)
xlabel('x (m)')
subplot(2,3,6)
xlabel('x (m)')

