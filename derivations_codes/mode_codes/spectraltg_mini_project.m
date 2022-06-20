%% TG equation for tutorial purposes.  Not tailored for instabilities
  clear all, close all
  format long, format compact
  N=200;
  [D,z]=cheb(N); D2=D^2; D2=D2(2:N,2:N); %differentiation matrics
  [zi w]=clencurt(N); % integration coefficients
  H=100; g=9.81; %basic dimensions
  dzpdzc=H*0.5; % Jacobian for derivatives in physical space
  bigl=H/2;  zphys=H-0.5*H*(z+1); %length scale and physical grid
  lambda=10000*H; % wave length of disturbance
  kphys=2*pi/lambda; % wave number of disturbance
  % combinations for near surface
  uamp=0.3; dj=0.1*H; % this pair of values has both mode-1 and mode-2 longwaves
  uamp=0.6; dj=0.1*H; % this pair of values has only mode-1 
  uamp=-0.6; dj=0.1*H; % this pair of values has only mode-1
  %combinations for near bottom
  uamp=0.3; dj=0.1*H; % this pair of values has both mode-1 and mode-2 longwaves
  uamp=1.25; dj=0.1*H; % this pair of values has only mode-1
  
  z0=H-0.25*H; d=0.05*H;
  dzd=1e-5*H;dzd2=dzd*dzd; %numerical differentiation parameters
  % background profiles
  my_density=@(z) 1-0.005*tanh((z-z0)/d);%
  my_d_density=@(z) (my_density(z+dzd)-my_density(z-dzd))/(2*dzd);;
  
  % dj sets the thickness near surface/bottom shear current
  %my_u=@(z) uamp*exp((z-H)/dj);
  % This is the version near the bottom
  my_u=@(z) uamp*exp(-z/dj);
  
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
E1pos=c1pos*phi1pos./(c1pos-uphysical);
[aa bb]=max(abs(E1pos));
E1pos=E1pos/E1pos(bb);
 
phi2pos=[0; Xim(1:N-1,sortindpos(2)); 0];
E2pos=c2pos*phi2pos./(c2pos-uphysical);
[aa bb]=max(abs(E2pos));
E2pos=E2pos/E2pos(bb);
 
[sorteeneg sortindneg]=sort(real(ee),'ascend');
c1neg=ee(sortindneg(1));
c2neg=ee(sortindneg(2));

phi1neg=[0; Xim(1:N-1,sortindneg(1)); 0];
E1neg=c1neg*phi1neg./(c1neg-uphysical);
[aa bb]=max(abs(E1neg));
E1neg=E1neg/E1neg(bb);
 
phi2neg=[0; Xim(1:N-1,sortindneg(2)); 0];
E2neg=c2neg*phi2neg./(c2neg-uphysical);
[aa bb]=max(abs(E2neg));
E2neg=E2neg/E2neg(bb);
 
figure(1)
clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
subplot(1,2,1)
plot(uphysical/max(abs(uphysical)),zphys,'k-',n2physical/n2max,zphys,'k--')
legend('U(z)','N^2(z)','Location','SouthEast')
hold on
plot([0 0],[-H 0],'k:')
plot([-2 2],-[z0 z0],'k:')
xlabel('Scaled U and N^2')
ylabel('z')
axis([-1.1 1.1 0 H])
grid on
title(['U_0 = ' num2str(uamp,4)])
subplot(1,2,2)
plot(E1pos,zphys,'k-',E1neg,zphys,'k--') % uncomment this if you want the
%more widely theoretically quoted E form of the eigenfunctions
%plot(phi1pos,zphys,'k-',phi1neg,zphys,'k--') % uncomment this if you want
%the eigenfunctions as given in the T-G equation

legend('positive','negative','Location','SouthEast')
hold on
plot([0 0],[0 H],'k:')
plot([-2 2],[z0 z0],'k:')
axis([-0.1 1.1 0 H])
grid on
xlabel('mode structure')
title([' (c1+,c1-)=(' num2str(c1pos,4) ',' num2str(c1neg,4) ')']);
figure(2)
clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
subplot(1,2,1)
plot(uphysical/max(abs(uphysical)),zphys,'k-',n2physical/n2max,zphys,'k--')
legend('U(z)','N^2(z)','Location','SouthEast')
hold on
plot([0 0],[0 H],'k:')
plot([-2 2],[z0 z0],'k:')
xlabel('Scaled U and N^2')
ylabel('z')
axis([-1.1 1.1 0 H])
grid on
title(['U_0 = ' num2str(uamp,4)])
subplot(1,2,2)
plot(E2pos,zphys,'k-',E2neg,zphys,'k--') % uncomment this if you want the
%more widely theoretically quoted E form of the eigenfunctions
%plot(phi2pos,zphys,'k-',phi2neg,zphys,'k--') % uncomment this if you want
%the eigenfunctions as given in the T-G equation
legend('positive','negative','Location','SouthEast')
hold on
plot([0 0],[0 H],'k:')
plot([-2 2],[z0 z0],'k:')
axis([-1.1 1.1 0 H])
grid on
xlabel('mode structure')
title([' (c2+,c2-)=(' num2str(c2pos,4) ',' num2str(c2neg,4) ')']);

