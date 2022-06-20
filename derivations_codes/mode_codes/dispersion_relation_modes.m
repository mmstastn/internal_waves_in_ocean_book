clear all, close all
% pseudospectral solution of the vertical structure of linear internal
% waves
%  lambda(-D^2+k^2 I)phi = N^2(z)phi
% on 0<z<H with phi(0)=phi(H)=0
  format long, format compact
  % This is what you need ot get the differentiation matrix and the grid
  N=200;
  [D,zc]=cheb(N); D2=D^2; D2=D2(2:N,2:N); 
   % here are the physical parameters and the scaling to the computational
   % domain [-1, 1]
  H=100; g=9.81;
  dzpdzc=H*0.5;
  dzcdzp=1/dzpdzc;
  % the physical parameters for the stratification
  a=0.005;
  z0=0.75*H;
  d=0.1*H;
  zphys=0.5*H*(zc+1);
  dz_num=1e-8*H;
  % inline functions for the density and the derivative of the density
  my_density=@(z) 1-a*tanh((z-z0)/d);
  my_d_density=@(z) (my_density(z+dz_num)-my_density(z-dz_num))/(2*dz_num);  
  n2physical=-g*my_d_density(zphys);
  n2max=max(n2physical);
  % wavelength and wavenumber
  cntr=0
  for pow=-1:0.125:8
      cntr=cntr+1;
      wavelength=H*10^pow; 
      kphys=2*pi/wavelength; kphys2=kphys^2;
      ks(cntr)=kphys;
      % make up the matrices for the e-val prog.

        %define B
        B=-D2*(1/dzpdzc)^2+kphys2*eye(size(D2));
        %define A
        A=diag(n2physical(2:end-1));
        % Solve the e-val prob 
        [ev ee]=eig(A,B);
        [cs csi]=sort(sqrt(diag(ee)),'descend');
        c1s(cntr)=cs(1);c2s(cntr)=cs(2);c3s(cntr)=cs(3);
        
        dk=1e-8*kphys;
        % make up the matrices for the e-val prog.
        %define B
        B=-D2*(1/dzpdzc)^2+(kphys+dk)^2*eye(size(D2));
        %define A
        A=diag(n2physical(2:end-1));
        % Solve the e-val prob 
        [ev ee]=eig(A,B);
        [cs csi]=sort(sqrt(diag(ee)),'descend');
        sig1p=cs(1)*(kphys+dk);sig2p=cs(2)*(kphys+dk);sig3p=cs(3)*(kphys+dk);
        % make up the matrices for the e-val prog.
        %define B
        B=-D2*(1/dzpdzc)^2+(kphys-dk)^2*eye(size(D2));
        %define A
        A=diag(n2physical(2:end-1));
        % Solve the e-val prob 
        [ev ee]=eig(A,B);
        [cs csi]=sort(sqrt(diag(ee)),'descend');
        sig1n=cs(1)*(kphys-dk);sig2n=cs(2)*(kphys-dk);sig3n=cs(3)*(kphys-dk);
        cg1s(cntr)=(sig1p-sig1n)/(2*dk);cg2s(cntr)=(sig2p-sig2n)/(2*dk);cg3s(cntr)=(sig3p-sig3n)/(2*dk);
  end
              figure(1), clf
        % This is for older versions of Matlab
%         set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
%                 'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
%                   'DefaultAxesFontWeight','bold');
        plot(ks*H,c1s/max(c1s),'kp-',ks*H,c2s/max(c1s),'ko--',ks*H,c3s/max(c1s),'k^:')
        hold on
        plot([min(ks) max(ks)]*H,[max(c2s) max(c2s)]/max(c1s),'k-')
        plot([min(ks) max(ks)]*H,[max(c3s) max(c3s)]/max(c1s),'k-')
        legend('mode-1','mode-2','mode-3','Location','Northeast')
        xlabel('k H')
        ylabel('phase speed, scaled')
        title('Phase Speed for the first three modes')
        axis([0 max(ks)*H 0 1])
        grid on
        
        figure(2), clf
        % This is for older versions of Matlab
%         set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
%                 'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
%                   'DefaultAxesFontWeight','bold');
        plot(ks*H,cg1s/max(c1s),'kp-',ks*H,cg2s/max(c1s),'ko--',ks*H,cg3s/max(c1s),'k^:')
        hold on
        plot([min(ks) max(ks)]*H,[max(c2s) max(c2s)]/max(c1s),'k-')
        plot([min(ks) max(ks)]*H,[max(c3s) max(c3s)]/max(c1s),'k-')
        legend('mode-1','mode-2','mode-3','Location','Northeast')
        xlabel('k H')
        ylabel('group speed, scaled')
        title('Group speed for the first three modes')
        axis([0 max(ks)*H 0 1])
        grid on