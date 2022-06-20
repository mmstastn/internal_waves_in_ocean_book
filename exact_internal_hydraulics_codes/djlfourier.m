%% This is the workhorse for solving the steady DJL by iteration
% This is the nonlinear piece based on the old eta
      forrh=dt*n2phys(zz-eta).*eta/(U0^2); forrhk=fft(forrh,[],1);
%Loop over Fourier modes
    for ki=1:N
% define the LH operator
      lhop=(1 + dt*k2(ki))*eye(Nz+1)-dt*Dzz;
      mytry=lhop*(etaf(ki,:)');max(abs(mytry));
% This imposes the BCs
      lhop(end,:)=0; lhop(end,end)=1;
      rh=etaf(ki,:)'+forrhk(ki,:)'; rh(end)=botbck(ki);
      lhop(1,:)=0; lhop(1,1)=1;
      %lhop(1,:)=etaf(ki,2)*(Dz(1,:)-1i*m*I(1,:))*exp(-1i*m*zz(ki,2)); 
      rh(1)=0;
      
      etaf(ki,:)=lhop\rh;   
      
    end
% end of loop over Fourier modes