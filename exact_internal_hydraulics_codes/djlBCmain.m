%% pseudospectral (cheb Fourier) iterative solver for supercritical DJL
close all
% The parameter FROM_SCRATCH needs to be set by user, default is 1
if exist('FROM_SCRATCH','var')==0
    FROM_SCRATCH=1
end
if FROM_SCRATCH
 clear all
 FROM_SCRATCH=1
 maxerr=1e5;
end
%% physical parameters
%U0=2.85; U02=U0*U0; % subcrit hill case
U0=3.0; U02=U0*U0; % supercrit hole case

% number of grid points (total is N+1) and the grid
L=2e3; H=100; h0=0.1*H; g=9.81;
% for a hole case
aa=0.02; z0=0.75*H; d=0.1*H; ishole=1;
% For a hill case
%aa=0.02; z0=0.25*H; d=0.1*H; ishole=0;
POL=-1+2*(1-ishole);
rhophys=@(s) (1-aa*tanh((s-z0)/d));
n2phys=@(s) (aa*g/d)*sech((s-z0)/d).^2;

% Grid points; start small and refine
N=192/2;
Nz=36/2; 
%% Main work horse loop
num_refine=1; % set the number of refinement cycles
for refi=0:num_refine
    N=N*2;Nz=Nz*2;
    make_grid
    % get the initial eta; I'm sure one could do better, but in practice this
    % seems OK
    if FROM_SCRATCH
     eta=zeros(N,Nz+1);
     eta0=POL*10*(sech(xx/2000).^4).*sech((zz-40)/40).*cos(0.5*pi*xx/L);
     
     eta=eta0;
     etaf=fft(eta,[],1);
     botbck=hk;
    end
    % set the quasi-time step (could be adjusted dynamically, but again this works)   
    dt=2e1;
    t=0;
    
    % Loop over BC updates
    for bcits=1:10
    % Loop over quasi-time steps  
     for jj=1:500
        t=t+dt;
        djlfourier
        etan=real(ifft(etaf,[],1)); etan=real(etan); diffnow=eta-etan;
        eta=etan;

       if (~mod(jj,100));
           djldiff1=zeros(size(eta));
           for dm=1:N
            djldiff1(dm,:)=Dzz*eta(dm,:)';
           end
           djldiff2=ifft(-k2mat.*etaf,[],1);
           djldiff3=n2phys(zz-eta).*eta/(U0^2);
           djldiff=djldiff1+djldiff2+djldiff3; djldiff=djldiff(:,2:end-1);
           djldiff=max(abs(djldiff(:)));
           disp(sprintf('%d %d %g %g',bcits,jj,max(abs(diffnow(:))),djldiff));
       end
     end
     % update the BCs
     botbcupdate
    % end of loop over bottom BCs
    end
    % sample plot
    figure(11);
    clf
    betterplots
    colormap((gray))
    subplot(2,1,1)
    contour(x,z,real(eta)',10);
    hold on
    plot(x(1,:),h,'k-','linewidth',2)
    hold off
    xlabel('x (m)')
    ylabel('z (m)')
    title('\eta (upper) and density (lower)')
    subplot(2,1,2)
    contour(x,z,rhophys(zz-eta)',10,'k')
    hold on
    plot(x(1,:),h,'k-','linewidth',2)
    hold off
    xlabel('x (m)')
    ylabel('z (m)')
    drawnow
    FROM_SCRATCH=0;
end

