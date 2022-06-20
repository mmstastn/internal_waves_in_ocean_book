%clear all, close all
  format long, format compact
  N=100;
  [D,z]=cheb(N); Dfull=D; D2=D^2; D2=D2(2:N,2:N); %differentiation matrics
  [zi w]=clencurt(N); % integration coefficients
  H=0.2; g=9.81; %basic dimensions
  dzpdzc=H*0.5; % Jacobian for derivatives in physical space
  bigl=H/2;  zphys=H-0.5*H*(z+1); %length scale and physical grid
  lambda=1000000*H; % wave length of disturbance
  kphys=2*pi/lambda; % wave number of disturbance
  uamp=0.0; % this value has both mode-1 and mode-2 longwaves
  %uamp=0.8; % this value has only mode-1 and the positive mode-2 longwave
  %uamp=-0.8; % this value has only mode-1 and the positive mode-2 longwave
  a_d=0.01; z0=H-0.05; d=0.02; 
  dzd=1e-5*H;dzd2=dzd*dzd; %numerical differentiation parameters
  % background profiles
  my_density=@(z) 1-a_d*tanh((z-z0)/d);%
  my_d_density=@(z) (my_density(z+dzd)-my_density(z-dzd))/(2*dzd);;
  
  my_u=@(z) uamp*z/H;
  my_uz=@(z) (my_u(z+dzd)-my_u(z-dzd))/(2*dzd); 
  my_uzz=@(z) (my_u(z+dzd)-2*my_u(z)+my_u(z-dzd))/dzd2;
  
  %create background profiles on grid
   rhophysical=my_density(zphys);
   n2physical=-g*my_d_density(zphys);
   uphysical=my_u(zphys);
   uzphysical=my_uz(zphys);
   uzzphysical=my_uzz(zphys);
   uzphysical(1)=uzphysical(2);
   uzzphysical(1)=uzzphysical(2);
   
  % some key parameters including Ri
  uzzmax=max(abs(uzzphysical));uzmax=max(abs(uzphysical));
  ri=n2physical./(uzphysical.*uzphysical); minri=min(ri);
  % time and velocity scales
  n2max=max(n2physical);  n0=sqrt(n2max); bigt=1/n0; bigu=bigl/bigt;
  
  % dimensionless profiles
  nn=n2physical(2:N)/n2max;
  u=uphysical(2:N)/bigu;
  uzz=uzzphysical(2:N)*bigl*bigl/bigu;
  k=kphys*bigl;
  k2=k*k;
  
  % make up the matrices for the e-val prog.
% eye(N) makes an N by N identity
% diag(vector) makes a diagonal matrix with the vector along the diagonal
p1=D2-k2*eye(N-1);
a11=diag(u)*p1-diag(uzz);
a12=-eye(N-1);
a21=diag(nn);
a22=diag(u);
b11=p1;
b22=eye(N-1);
%define A
A=zeros(2*(N-1),2*(N-1));
A(1:N-1,1:N-1)=a11;
A(1:N-1,N:2*(N-1))=a12;
A(N:2*(N-1),1:N-1)=a21;
A(N:2*(N-1),N:2*(N-1))=a22;
%define B
B=zeros(2*(N-1),2*(N-1));
B(1:N-1,1:N-1)=b11;
B(N:2*(N-1),N:2*(N-1))=b22;

% Solve the e-val problem
% and this is just plain eig
[Xim eedim]=eig(A,B);

% Post-process the eigenvalues when there are no imaginary parts

ee=diag(eedim)*bigu; %redimensionalize
[sorteepos sortindpos]=sort(real(ee),'descend');
c1pos=ee(sortindpos(1));
c2pos=ee(sortindpos(2));

phi1pos=[0; Xim(1:N-1,sortindpos(1)); 0];

 
phi1_z=Dfull*phi1pos;
mymax=max(phi1_z);
for ii=1:4
  uwnl(:,ii)=phi1_z*max(uisw(:,ii))/mymax;
end


figure(1)
clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
      hold on
      
      for ii=1:4
          subplot(2,2,ii)
          plot(uisw(:,ii),zisw,'k')
          hold on
          h(ii)=plot(uwnl(:,ii),zphys-H,'--');
          set(h(ii),'Color',[0.6 0.6 0.6])
          legend('DJL','WNL','Location','SouthEast')
          grid on
          axis([-0.025 0.04 -0.2 0])
      end
      subplot(2,2,1)
      ylabel('z (m)')
      subplot(2,2,3)
      ylabel('z (m)')
      xlabel('u (m/s)')
      subplot(2,2,4)
      xlabel('u (m/s)')
      
