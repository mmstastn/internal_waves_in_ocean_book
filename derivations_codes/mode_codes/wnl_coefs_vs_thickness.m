clear all, close all
% pseudospectral solution of the vertical structure of linear internal
% waves
%  lambda(-D^2+k^2 I)phi = N^2(z)phi
% on 0<z<H with phi(0)=phi(H)=0
  format long, format compact
  % This is what you need ot get the differentiation matrix and the grid
  N=200;
  [D,zc]=cheb(N); D2=D^2; D2=D2(2:N,2:N); 
  [zi w]=clencurt(N);
   % here are the physical parameters and the scaling to the computational
   % domain [-1, 1]
  H=100; g=9.81;
  dzpdzc=H*0.5;
  dzcdzp=1/dzpdzc;
  % the physical parameters for the stratification
  a=0.005;
  z0=0.75*H;
  zphys=0.5*H*(zc+1);
  dz_num=1e-8*H;
  cntr=0;
  for d0frac=0.01:0.01:0.99
      cntr=cntr+1;
      d=H*d0frac;
  % inline functions for the density and the derivative of the density
      my_density=@(z) 1-a*tanh((z-z0)/d);
      my_d_density=@(z) (my_density(z+dz_num)-my_density(z-dz_num))/(2*dz_num);  
      n2physical=-g*my_d_density(zphys);
      n2max=max(n2physical);
  
      % make up the matrices for the e-val prog.

        %define B
        B=-D2*(1/dzpdzc)^2;
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
        phi1p=D*[0;phi1;0]*dzcdzp;
        S1=sum(w'.*(phi1p.^2)*dzpdzc);
        r10_1(cntr)=-0.75*sum(w'.*(phi1p.^3)*dzpdzc)/S1;
        r01_1(cntr)=-0.5*c1*sum(w'.*([0;phi1;0].^2)*dzpdzc)/S1;
        phi2p=D*[0;phi2;0]*dzcdzp;
        S2=sum(w'.*(phi2p.^2)*dzpdzc);
        r10_2(cntr)=-0.75*sum(w'.*(phi2p.^3)*dzpdzc)/S1;
        r01_2(cntr)=-0.5*c2*sum(w'.*([0;phi2;0].^2)*dzpdzc)/S2;
        c1s(cntr)=c1;c2s(cntr)=c2;
        d0s(cntr)=d;
  end
  figure(1), clf
  set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
            'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
              'DefaultAxesFontWeight','bold');
  refnum=15        
  plot(c1s/c1s(refnum),d0s/H,'bo',r10_1/r10_1(refnum),d0s/H,'ko',r01_1/r01_1(refnum),d0s/H,'ro')
  hold on
  plot([-2.5 2.5],[1 1]*d0s(refnum)/H)
  ylabel('pycnocline thickness scaled by H')
  xlabel('scaled c, r_{10}, r_{01} for mode 1')
  legend('c','r_{10}','r_{01}','Location','NorthEast')
  axis([-0.5 2 0 1])
  grid on
  figure(2), clf
  set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
            'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
              'DefaultAxesFontWeight','bold');
  plot(c2s/c2s(refnum),d0s/H,'bo',r10_2/r10_2(refnum),d0s/H,'ko',r01_2/r01_2(refnum),d0s/H,'ro')
  hold on
  plot([-2.5 2.5],[1 1]*d0s(refnum)/H)
  axis([-0.5 2 0 1])
  ylabel('pycnocline thickness scaled by H')
  xlabel('scaled c, r_{10}, r_{01} for mode 2')
  legend('c','r_{10}','r_{01}','Location','NorthEast')
  grid on
  
      