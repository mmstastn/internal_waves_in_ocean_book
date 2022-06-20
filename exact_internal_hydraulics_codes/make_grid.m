if refi>0
    zzold=zz;xxold=xx;
end
[D0,z0]=cheb(Nz); D02=D0*D0;
if (ishole)
 Hz=H+h0; Hz2=Hz*Hz; z=0.5*(1+z0)*Hz-h0; Dz=(2/Hz)*D0; Dzz=(4/Hz2)*D02; 
else
  Hz=H; Hz2=Hz*Hz; z=0.5*(1+z0)*Hz; Dz=(2/Hz)*D0; Dzz=(4/Hz2)*D02; 
end
x0=linspace(-1,1,N+1);
x0=x0(1:end-1);
x=x0*L; dx=x(2)-x(1); 
[zz xx]=meshgrid(z,x);
dk=pi/L;
ksvec=zeros(size(x));
ksvec(1)=0; ksvec(N/2+1)=0;
for ii=2:(N/2)
   ksvec(ii)=ii-1;
   ksvec(N/2+ii)=-N/2 + ii -1;
end
k=ksvec*dk; k2=k.*k; ik=sqrt(-1)*k;
k2mat=repmat(k2',1,Nz+1);
if (ishole)
 %h=h0*(-(sech((x-10e2)/50)-0.25*sech(x/200)+sech((x+10e2)/50)).*cos(0.5*pi*x/L));
 h=-h0*sech(x/450).^2;
 %h=-h0*sech((x+500)/450).^2-h0*sech((x-500)/450).^2
else
    %h=h0*sech((x+500)/50).^2+h0*sech((x-500)/50).^2;
    h=h0*sech(x/450).^2;
end
 % a=5e3;
% h=h0*a^2./(x.^2+a^2);
hk=fft(h);
if refi>0
    etanew=interp2(xxold',zzold',eta',xx',zz','spline');
    botbcold=real(ifft(botbck));
    botbc=interp1(xxold(:,1),botbcold,xx(:,1),'spline')';
    botbck=fft(botbc);
    figure(13)
    clf
    betterplots
    subplot(2,1,1)
    plot(xxold(:,1),botbcold,'b-',xx(:,1),botbc,'r-')
    subplot(2,2,3)
    pcolor(xxold,zzold,eta),shading flat
    subplot(2,2,4)
    pcolor(xx,zz,etanew'),shading flat
    drawnow
    eta=etanew';etaf=fft(eta,[],1);
end
