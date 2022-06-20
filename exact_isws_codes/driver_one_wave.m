%% One lab scale wave using our DJL solver
clear all,close all
verbose=0; %If this is set to one you will see all the progress during solving
picshift=0; % This is used in the picture maker for multiple cases
A  = 4e-5; % APE for wave
NX = 256; % grid 
NZ = 256; % grid
H  = 0.2; %depth
L  = 2; %width

start_time = clock;
dzd=H*1e-4;dzd2=dzd*dzd;
uamp=0.0;
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
iswpic_one_wave_bw

end_time=clock;
fprintf('Total wall clock time: %f\n',etime(end_time, start_time));
