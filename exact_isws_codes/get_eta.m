%%%%%% this is the main iteration loop get_eta

% initialise t for recording timing data
t.init = 0; t.solver = 0; t.plot = 0; t.loop = 0; t.int1 = 0; tic

% set min and max number of iterations, g
min_iterate=25;
max_iterate=1000;
g=9.81;

% number of digits to solve to, 1e-5 means 5 
epsilon = 1e-4;
% number of steps for Simpson's rule integration
NS=20;


%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare grid, etc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx= (L*(2*(1:NX)-1) / (2*NX));
zz=-(H*(2*(1:NZ)-1) / (2*NZ));
%zz=-(H*(2*(NZ:-1:1)-1) / (2*NZ));

[x, z] = meshgrid(xx,zz);
dx=abs(xx(2)-xx(1));
dz=abs(zz(2)-zz(1));
dxdz=dx*dz;

ks = (pi/L) * [0:(NX-1) -NX:-1];
% this minus is to account for the ordering of the zz s
ms = -(pi/H) * [0:(NZ-1) -NZ:-1]';
%ms = (pi/H) * [0:(NZ-1) -NZ:-1]';

KS = repmat(ks, [2*NZ 1]);
MS = repmat(ms, [1 2*NX]);
LAP = -KS.^2  -MS.^2;
ONE_OVER_LAP = 1./LAP;
ONE_OVER_LAP(LAP == 0)=0;

clear ms ks

%%%%% Prepare initial guess and iteration independant steps %%%%%%%%%%%%%%%

% this is loop-independant
den_fine_z = md_density(z);

if ( exist('eta0', 'var') )
    % Initial guess provided, interpolate if needed
    if( ~isequal(x,x0) && ~isequal(z,z0))
        % use FFT interpolation
        [NZ0 NX0] = size(eta0);
        RX=NX/NX0; RZ=NZ/NZ0;

        eta0in = [eta0 -fliplr(eta0)];
        eta0out = interpft(eta0in, 4*NX, 2);
        eta0out = circshift(eta0out, [0 (RX-1)]);
        eta0 = eta0out(:, 1:2:2*NX);

        eta0in = [eta0; -flipud(eta0)];
        eta0out = interpft(eta0in, 4*NZ, 1);
        eta0out = circshift(eta0out, [(RZ-1) 0]);
        eta0 = eta0out(1:2:2*NZ, :);
        
        clear x0 z0 eta0in eta0out
        
        % increasing resolution by FFT interpolation
        % usually gives an extremely good guess, so
        % reduce min_iterate to save time
        min_iterate = 1; 
    end
else
    % Initial guess by weakly nonlinear theory

    Dz = md_diff(dz, NZ, 1, 'not periodic');
    Dzz = md_diff(dz, NZ, 2, 'not periodic');
    Dzzc = Dzz(2:end-1, 2:end-1);

    % get n2, u and uzz data
    tmpz   =  zz(2:end-1)';
    n2vec  = -g*md_d_density(tmpz);
    uvec   =  md_u(tmpz);
    uzzvec =  md_uzz(tmpz);

    % create diagonal matrices
    N2  = diag(n2vec );
    U   = diag(uvec  );
    Uzz = diag(uzzvec);

    % setup quadratic eigenvalue problem
    A0 = sparse(N2 + U.*U*Dzzc - U*Uzz);
    A1 = sparse(-2*U*Dzzc + Uzz);
    A2 = sparse(Dzzc);

    % solve eigenvalue problem; extract first eigenvalue&eigenmode
    [V,cc]=polyeig(A0,A1,A2);
    [c ii]=sort(cc,'descend');
    c0=c(1);

    % add boundary conditions
    phi1=[0; V(:,ii(1)); 0];
    uvec=[0; uvec; 0];
    
    % compute E1, normalise
    E1=c0*phi1./(c0-uvec);
    E1=E1/max(abs(E1));
    E1=abs(E1);

    % compute r10 and r01
    E1p=Dz*E1; E1p2=E1p.^2; E1p3=E1p.^3;
    bot=sum((c0-uvec).*E1p2);
    r10=(-0.75/c0)*sum((c0-uvec).*(c0-uvec).*E1p3)/bot;
    r01=-0.5*sum((c0-uvec).*(c0-uvec).*E1.*E1)/bot;
if(verbose)
    fprintf('WNL gives: c_lw = %f, r10 = %f, r01=%f\n\n',c0, r10, r01);
end
    % now optimise the b0, lambda parameters
    E = repmat(E1, 1, NX); E(:,1)=0; E(:,end)=0;
    lambda=L/8;
    ampl=sign(r10)*0.15*H;
    V = (1+(2/3)*r10*ampl)*c0;
    if(verbose)
    fprintf('init ampl = %f, lambda = %f, V=%f\n',ampl, lambda, V);
    end
    flag = 1;
    while(flag)
      B = ampl*sech((x-L/2)/lambda).^2;
      etai = B.*E;
      etai(1,:)=0; etai(:,1)=0; etai(:,end)=0; etai(end,:)=0;

%       figure(1)
%       subplot(2,2,1); pcolor(x,z,psi); colorbar; title('Initguess'); pause(0.01);
%       drawnow;
%       subplot(2,2,2); plot(E1,zz)

      % use Simpson's rule here to find F
      den_value = md_density(z-etai);
      ff        = den_fine_z + den_value;
      for n=(1:2:NS-1)/NS
          ff = ff + 4*md_density(z-n*etai);
      end
      for n=(2:2:NS-2)/NS
          ff = ff + 2*md_density(z-n*etai);
      end
     
      ff        = (-ff/(3*NS) + den_value).*etai;
      F  = dxdz*sum((ff(:))) /H;
     
      % new ampl by rescaling
      afact = A / F;
      afact = min(afact, 1.5);
      afact = max(afact, 0.75);
      ampl = ampl * afact;

      % new lambda
      lambda = sqrt (  -6*r01 / (c0 * r10 * ampl)  );
      if ~isreal(lambda)
          disp('problem finding new lambda !!')
      end
      
      % new V
      V = (1+(2/3)*r10*ampl)*c0;
if(verbose)
      fprintf('F=%e, desired = %e, rescaling ampl by factor of %f...\n',F, A,afact);
      fprintf('new  ampl = %f, lambda = %f, V=%f\n\n',ampl, lambda, V);
end
      % we're close enough now or an initial guess
% Marek     if abs(afact-1) < 0.01
          eta0=etai;
          flag=0;
 % Marek    end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!! Something weird happened, so negate the init guess %
%eta0=-eta0;                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Dz Dzz Dzzc phi1 psi N2 U Uzz A0 A1 A2 B E 
end


% check if there is a nonzero background profile, save time if it's zero
uflag = any(md_u(zz));

% prepare the first iteration
maxeta     = max((eta0(:)));
mineta     = min((eta0(:)));
lambda0    = g*H/(c0*c0);
if(verbose)
fprintf('Initial guess: maxeta = %1.6e,   mineta = %1.6e\n',maxeta,mineta);
end
flag = 1; iterate_num = 0; t.init = toc; tic; % record timing data

while(flag)
    iterate_num = iterate_num + 1; 
    
    % compute S
    S = -md_d_density(z-eta0).*eta0/H;
    
    % compute R, assemble RHS
    if(uflag)
        array = [eta0 -fliplr(eta0); -flipud(eta0) fliplr(flipud(eta0))];    
        ARRAY = fft2(array);
        array = ifft2(i*KS.*ARRAY); eta0x = array(1:NZ, 1:NX);
        array = ifft2(i*MS.*ARRAY); eta0z = array(1:NZ, 1:NX);
        uhat  = md_u (z-eta0)/c0;
        uhatz = md_uz(z-eta0)/c0;
        R = (uhatz ./ (uhat -1)) .* ( 1 - (eta0x.^2 + (1-eta0z).^2));
        rhs = - ( lambda0 * (S./((uhat-1).^2))  + R );
    else
        rhs = - lambda0 * S;
    end
    
    % record timing data
    tmp = toc; t.loop = t.loop + tmp; t.xloop = tmp; tic;

    % solve the linear Poisson problem
    array = [rhs -fliplr(rhs); -flipud(rhs) fliplr(flipud(rhs))];
    array = ifft2(fft2(array) .* ONE_OVER_LAP);
    nu = array(1:NZ, 1:NX);

    % record timing data
    tmp = toc; t.solver = t.solver + tmp; t.xsolver = tmp; tic;

    % use Simpson's rule here to find ff
    den_value = md_density(z-eta0);
    ff        = den_fine_z + den_value;
    for n=(1:2:NS-1)/NS
        ff = ff + 4*md_density(z-n*eta0);
    end
    for n=(2:2:NS-2)/NS
        ff = ff + 2*md_density(z-n*eta0);
    end
    ff        = (-ff/(3*NS) + den_value).*eta0;

    % record timing data
    tmp = toc; t.int1 = t.int1 + tmp; t.xint1 = tmp; tic;

    % compute F, S1, S2
    F  = dxdz * sum( ff(:)         )/H;
    S1 = dxdz * sum( S(:).*nu(:)   );
    S2 = dxdz * sum( S(:).*eta0(:) );

    % new lambda
    lambda = real(lambda0*(A-F+S2)/S1); 

    % check if lambda is OK
    if(lambda < 0)
        if(verbose)
        fprintf('new lambda = %1.6e\n',lambda); 
        fprintf('new lambda has wrong sign ==> nonconvergence of iteration procedure.\n'); 
        fprintf('   A = %1.8e, F = %1.8e\n',A,F); 
        fprintf('   S1 = %1.6e, S2 = %1.6e, S2/S1 = %1.6e\n',S1,S2,S2/S1);  
        end
        break;
    end

    % compute new c, eta, etax, etaz
    c    = sqrt(g*H/lambda);
    % Marek: I have added some underrelaxation here for low Ri runs
    % This is what Kevin used for his low Ri waves
    eta  = 0.8*eta0+0.2*(lambda/lambda0) * real(nu);

    % compute maximum difference from previous eta
    max_diff = max(abs(eta(:)-eta0(:)));
    maxeta = max((eta(:)));
    mineta = min((eta(:)));
    
    % compute residuals
%     residual = md_residual(LAP, eta, etax, etaz, z, g, c);
%     max_residual=max(abs(residual(:)));
    
    % report on state of the operation
    if(verbose)
    fprintf('\niterate_num = %4d:  max_diff = %1.8e\n',iterate_num,max_diff); 
%    fprintf('                 max_residual = %1.8e\n',max_residual); 
    fprintf('      maxeta  = %1.8e,  mineta = %1.8e\n',maxeta,mineta); 
    fprintf('      lambda  = %1.6e,     c  = %1.6e\n',lambda,c); 
    fprintf('      A  = %1.6e,  F = %1.6e\n',A,F); 
    fprintf('      S1 = %1.6e, S2 = %1.6e, S2/S1 = %1.6e\n',S1,S2,S2/S1);  
    end
    % stop conditions
    if( iterate_num >= min_iterate)
        if ( abs(max_diff / max([abs(maxeta),abs(mineta)])) < epsilon)
            flag = 0;
        end
        if(iterate_num >= max_iterate)
            flag = 0; 
            fprintf('******* Reached maximum number of iterations *******\n'); 
        end
    end

    % record timing data
    tmp = toc; t.loop = t.loop + tmp; t.xloop = t.xloop + tmp; tic;

%     % show progress
%     if(iterate_num < 2)% | (0 == mod(iterate_num, 3))
%         % compute residuals for plot
%         residual = md_residual(LAP, KS, MS, eta, z, g, c);
%         md_showprogress(eta, rhs, residual, x, z, iterate_num, A);
%     end
 
    % iteration shift
    lambda0  = lambda;
    c0       = c;
    eta0     = eta;
   
    % record timing data
    tmp = toc; t.plot = t.plot + tmp; t.xplot = tmp; tic;
 if(verbose)
    fprintf('Solver time:   %f seconds\n', t.xsolver);
    fprintf('Integral time: %f seconds\n', t.xint1);
    fprintf('Loop time:     %f seconds\n', t.xloop);
    fprintf('Plot time:     %f seconds\n', t.xplot);
    fprintf('----------------------------------------------------------\n');
 end
end  %% end iteration loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = toc; t.loop = t.loop + tmp; t.xloop = t.xloop + tmp; tic;

residual=0;
tmp = toc; t.plot = t.plot + tmp; t.xplot = tmp;


t.total = t.init + t.loop + t.int1 + t.solver + t.plot;
if(verbose)
fprintf('Initialisation time:  %7.2f seconds\n', t.init);
fprintf('Total Loop time:      %7.2f seconds\n', t.loop);
fprintf('Total Integral time:  %7.2f seconds\n', t.int1);
fprintf('Total Solver time:    %7.2f seconds\n', t.solver);
fprintf('Total Plotting time:  %7.2f seconds\n', t.plot);
fprintf('Total Time:           %7.2f seconds\n', t.total);
end
%whos
vars=whos; 
memtotal=0;
for ii=1:length(vars)
  memtotal = memtotal +vars(ii).bytes;
end
%  cfprintf('Memory usage: %.2f MB (%.2fx sizeof(eta))\n', ...
%           memtotal/1024.^2,memtotal/(8*numel(eta)));

fprintf('Finished [NZ,NX]=[%dx%d], A=%1.6e \n', NX,NZ,A);