% Create the temporary wavenumbers
ks = (pi/L) * [0:(NX-1) -NX:-1];
% the minus is to account for the ordering of the z grid
ms = -(pi/H) * [0:(NZ-1) -NZ:-1]';

KS = repmat(ks, [2*NZ 1]);
MS = repmat(ms, [1 2*NX]);
       
% extend the variable of interest
array = [eta -fliplr(eta); -flipud(eta) fliplr(flipud(eta))];    
ARRAY = fft2(array);
array = ifft2(i*KS.*ARRAY); etax = array(1:NZ, 1:NX);
array = ifft2(i*MS.*ARRAY); etaz = array(1:NZ, 1:NX);

up  = md_u (z-eta);
uu  = md_u (z);

%velocities
u=real(up+(c-up).*etaz);
uwave=u-c;

w=real(up.*(etax)-c*etax);

%kinteic energy
kewave=0.5*(uwave.^2 + w.^2);

%vorticity
array = [kewave fliplr(kewave); flipud(kewave) fliplr(flipud(kewave))];    
ARRAY = fft2(array);
array = ifft2(i*KS.*ARRAY); kewavex = real(array(1:NZ, 1:NX));
array = ifft2(i*MS.*ARRAY); kewavez = real(array(1:NZ, 1:NX));

array = [u fliplr(u); flipud(u) fliplr(flipud(u))];    
ARRAY = fft2(array);
array = ifft2(i*KS.*ARRAY); ux = real(array(1:NZ, 1:NX));
array = ifft2(i*MS.*ARRAY); uz = real(array(1:NZ, 1:NX));

array = [w -fliplr(w); -flipud(w) fliplr(flipud(w))];    
ARRAY = fft2(array);
array = ifft2(i*KS.*ARRAY); wx = real(array(1:NZ, 1:NX));
array = ifft2(i*MS.*ARRAY); wz = real(array(1:NZ, 1:NX));

vort=uz-wx;

%density and gradient Richardson number
den=md_density(z-eta);
n2p=-9.81*md_d_density(z-eta);
ri=n2p./(uz.*uz);
rifilt=ri.*(n2p>1e-4&abs(u)>1e-3*c)+100.*(n2p<1e-4|abs(u)<1e-3*c);
uv=u(:,NX/2);
uzv=uz(:,NX/2);

riv=ri(:,NX/2);
n2v=n2p(:,NX/2);
zv=z(:,NX/2);

% Filter the Richardson number since having Ri<0.25 is irrelevant if there
% is no stratification (N^2 is too small)
rifilt=ri+100*(n2p<1e-6);
rivfilt=riv+100*(n2v<1e-6);

% The second invariant which contains shear
i2=-wz.^2-0.25*(uz+wx).^2;
i2absmax=max(abs(i2(:)));

figure(2)