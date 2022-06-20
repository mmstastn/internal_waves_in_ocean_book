clear all,close all
verbose=0;
picshift=0;

A  = 2e-2; % APE for wave sets the amplitude of the wave implicitly
NX = 256; % grid 
NZ = 256; % grid
H  = 20.0; %depth
L  = 1000.0; %width

start_time = clock;
dzd=H*1e-4;dzd2=dzd*dzd;
uamp=0.0;
from_data=1;
% This part uses data for the stratification
if from_data==1
 load sample_data.csv
 zd=sample_data(:,1)-20;
 dd=sample_data(:,3);
 md_density=@(z) interp1(zd,dd,z,'spline')/sample_data(1,3);
 md_d_density=@(z) (md_density(z+dzd)-md_density(z-dzd))/(2*dzd);
else
    % This part specifies the stratification analytically
 a_d=0.02; z0_d=0.15; d_d=0.005;  % for wave of elevation
% a_d=0.02; z0_d=0.05; d_d=0.005;  % for wave of depression
 md_density=@(z) 1-a_d*tanh((z+z0_d)/d_d);
 md_d_density=@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;
end 
% This is the background current if there is any
zj=0.5*H;dj=0.4*H;
md_u=@(z) uamp*tanh((z+zj)/dj);
md_uz=@(z) (md_u(z+dzd)-md_u(z-dzd))/(2*dzd);
md_uzz=@(z) (md_u(z+dzd)-2*md_u(z)+md_u(z-dzd))/dzd2;
% solve problem once
get_eta
% plot things
iswpost
iswpic

% %increase resolution
 x0=x; z0=z;
 NX=NX*2; NZ=NZ*2;
 get_eta
% get the various secondary fields like velocities, vorticity, density, N^2
% and Ri
iswpost
iswpic

end_time=clock;
fprintf('Total wall clock time: %f\n',etime(end_time, start_time));
