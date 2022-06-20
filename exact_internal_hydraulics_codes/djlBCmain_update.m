%% pseudospectral (cheb Fourier) iterative solver for supercritical DJL
%% loop over U0
clear all,close all
numUs=6; % the number of U values
plotrows=2,plotcols=3; % make sure these multiply to give numUs
Us=linspace(3.0,3.5,numUs);
maxerr=1e5;
% initialize graphics
figure(11);
clf
betterplots
colormap((gray))
for cntr=1:numUs
    FROM_SCRATCH=1
    maxerr=1e5
      % This uses present solution as a guess
    if cntr==1
        FROM_SCRATCH=1;
    else
        FROM_SCRATCH=0;
    end
    %% physical parameters
    %U0=10; N0=0.01;N02=N0*N0; U02=U0*U0;
    U0=Us(cntr); U02=U0*U0;
    % number of grid points (total is N+1) and the grid
    L=2e3; H=100; h0=0.1*H; g=9.81;
    aa=0.02; z0=0.75*H; d=0.1*H; ishole=1;
    POL=-1+2*(1-ishole);
    rhophys=@(s) (1-aa*tanh((s-z0)/d));
    n2phys=@(s) (aa*g/d)*sech((s-z0)/d).^2;

    % Grid points; start small and refine
    N=256;
    Nz=36; 
    %% Main work horse loop
    
    refi=0
    make_grid
    

    % get the initial eta; I'm sure I could do better, but in practice this
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
    
    mxabseta(cntr)=max(abs(eta(:)));
    % sample plot
    figure(11)
    subplot(plotrows,plotcols,cntr)
    contour(x,z,rhophys(zz-eta)',10,'k')
    hold on
    plot(x(1,:),h,'k-','linewidth',2)
    hold off
    xlabel('x (m)')
    ylabel('z (m)')
    drawnow
    
end
