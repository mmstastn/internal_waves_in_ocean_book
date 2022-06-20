%% This is version of the DJL code starts with a data set
%% (the data is from Monterrey Bay; courtesy of Ryan Walter)
%% This specifies background density and currents,
%% but the waves can propagate in different directions
%% and hence feel different background currents.
clear all,close all

verbose=0;

A  = 1.5*1.2e-2; % APE for wave sets the amplitude of the wave implicitly
NX = 512; % grid 
NZ = 128; % grid
H  = 20; %depth
L  = 1000; %width

start_time = clock;
dzd=H*1e-4;dzd2=dzd*dzd;
uamp=0.05;
from_data=1;
% This part uses data for the stratification
if from_data==1
 load upstream_profiles
 zd=zv-20;
 dd=myrho1;
 uprofile=myu1;%uu(1:15)=uu(16);
 vprofile=myv1; %vv(1:11)=vv(12);
 
 md_density=@(z) interp1(zd,dd,z,'spline')/dd(1);
 md_d_density=@(z) (md_density(z+dzd)-md_density(z-dzd))/(2*dzd);
 set_u=@(z) 0.5*interp1(zd,uprofile,z,'spline')
 set_v=@(z) 0.5*interp1(zd,vprofile,z,'spline')
else
    % This part specifies the stratification analytically
 a_d=0.02; z0_d=0.05; d_d=0.01;  % for wave of elevation
% a_d=0.02; z0_d=0.05; d_d=0.005;  % for wave of depression
 md_density=@(z) 1-a_d*tanh((z+z0_d)/d_d);
 md_d_density=@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;

% This is the background current if there is any
zj=0.5*H;dj=0.1*H;
set_u=@(z) uamp*tanh((z+zj)/dj);
set_v=@(z) 0*tanh((z+zj)/dj);
end 
% Which angles do you want to plot
numangs=4;
min_ang=pi/16;max_ang=pi/2; % for showing the one odd case
%min_ang=0;max_ang=pi; % for showing with and against


angs=linspace(min_ang,max_ang,numangs);
figure(1)
clf
colormap darkjet
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',14,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',14,...
      'DefaultAxesFontWeight','bold');
figure(2)
clf
colormap darkjet
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',14,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',14,...
      'DefaultAxesFontWeight','bold');
  figure(3)
clf
colormap darkjet
set(gcf,'DefaultLineLineWidth',3,'DefaultTextFontSize',14,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',14,...
      'DefaultAxesFontWeight','bold');
  % This is the loop over the angles
for angsi=1:numangs
    ang=angs(angsi)
    md_u=@(z) cos(ang)*set_u(z)+sin(ang)*set_v(z);
    md_uz=@(z) (md_u(z+dzd)-md_u(z-dzd))/(2*dzd);
    md_uzz=@(z) (md_u(z+dzd)-2*md_u(z)+md_u(z-dzd))/dzd2;
    % solve problem once
    get_eta
    if(iterate_num >= max_iterate)
            fprintf('******* Reached maximum number of iterations *******\n'); 
    end
    % plot things
    iswpost
    cs(angsi)=c;
    figure(1)
    % If you change the number of angles you will want to change
    % the subplot information
    hnow=subplot(2,2,angsi);
    [conts,h]=contourf(x,z,den,20);
    set(h,'LineColor','none')
    colormap darkjet
    figure(2)
    subplot(2,2,angsi)
    plot(u(:,NX/2),z(:,NX/2),'b-',u(:,1),z(:,1),'r-')  
    axis([-0.2 0.2 -H 0])
    grid on
    figure(3)
    subplot(2,2,angsi)
    plot(ri(:,end/2),z(:,end/2),'b-',ri(:,1),z(:,1),'r-')
    axis([-0.1 0.9 -H 0])
    grid on
    x0=x; z0=z;
end
figure(2)
subplot(2,2,1)
title('u wave-induced blue')
subplot(2,2,2)
title('u background red')
figure(3)
subplot(2,2,1)
title('Ri wave-induced blue')
subplot(2,2,2)
title('Ri background red')

end_time=clock;
fprintf('Total wall clock time: %f\n',etime(end_time, start_time));
