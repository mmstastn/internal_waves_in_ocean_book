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
  for pow=-2:0.1:2
      wavelength=H*10^pow; 
      kphys=2*pi/wavelength; kphys2=kphys^2;

      % make up the matrices for the e-val prog.

        %define B
        B=-D2*(1/dzpdzc)^2+kphys2*eye(size(D2));
        %define A
        A=diag(n2physical(2:end-1));
        % Solve the e-val prob 
        [ev ee]=eig(A,B);
        [cs csi]=sort(sqrt(diag(ee)),'descend');
        c1=cs(1);c2=cs(2);c3=cs(3);
       % This makes sure that the eigenfunction has a maximum of 1
        phi1=ev(:,csi(1)); 
        mxphi1=max(phi1);
        mnphi1=min(phi1);
        mxabs=max(abs(phi1));
            if abs(mnphi1)==mxabs
                phi1=-phi1/mxabs;
            else
                phi1=phi1/mxabs;
            end
        phi2=ev(:,csi(2)); 
        mxphi2=max(phi2);
        mnphi2=min(phi2);
        mxabs2=max(abs(phi2));
            if abs(mnphi2)==mxabs2
                phi2=-phi2/mxabs2;
            else
                phi2=phi2/mxabs2;
            end
        figure(1), clf
        % This is for older versions of Matlab
%         set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
%                 'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
%                   'DefaultAxesFontWeight','bold');
        plot(phi1,zphys(2:end-1),'b',phi2,zphys(2:end-1),'k',n2physical/n2max,zphys,'r' )
        legend('\phi^{(1)}','\phi^{(2)}','N^2','Location','Southeast')
        xlabel('\phi and N^2')
        ylabel('z')
        title(['Wavelength scaled by H =' num2str(wavelength/H,3) ', c = ' num2str(c1,3)])
        axis([-1 1 0 H])
        grid on
        pause
  end
      