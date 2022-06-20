%% This script updates the boundary condition for the steady DJL

% Interpolate to the actual boundary at each x_i
  for ki=1:N
     etab(ki)=interp1(z,eta(ki,:),h(ki),'spline');
%     etab(ki)=chebint(eta(ki,:),h(ki));
  end
% Compute the error
 en=etab-h;
 newmaxerr=max(abs(en));
 disp(['maxerr = ' num2str(maxerr) ' newmaxerr = ' num2str(newmaxerr)])
 enf=fft(en);
 figure(1)
 clf
 betterplots
 subplot(2,1,1)
 plot(x,etab,'b-',x,h,'r-')
 title('eta at z=h in blue, h in red')
 subplot(2,1,2)
 plot(k,abs(fft(etab)),'b.',k,abs(fft(h)),'r.')
 title('Spectra of eta at z=h in blue, h in red')
 drawnow
% define the new bottom BC
%if(newmaxerr<maxerr)
 disp('updating BCs')
 botbck=botbck-enf;
 maxerr=newmaxerr;
%end
