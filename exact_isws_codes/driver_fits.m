%% lab scale waves using our DJL solver for fitting
clear all,close all
verbose=0; %If this is set to one you will see all the progress during solving
As=[1 5 10 15]*1e-6;
figure(1)
clf
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
      'DefaultAxesFontWeight','bold');
colormap gray
for casei=1:length(As)
    uamp=0;
    A  = As(casei); % APE for wave
    NX = 512; % grid 
    NZ = 512; % grid
    H  = 0.2; %depth
    L  = 3.5; %width

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
    uisw(:,casei)=u(:,NX/2);
end
zisw=z(:,1);
uvert_fit