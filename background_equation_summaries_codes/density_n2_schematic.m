% Density and N^2 schematic
z=linspace(0,100,1001);
figure(1)
clf
betterplots
z=linspace(0,100,1001);
rhobarlab=1023*(1+0.005*tanh((z-80)/2.5));
rhobarocean=1023*(1+0.005*tanh((z-80)/5)+0.0005*tanh((z-40)/30));
n2lab=(0.005/2.5)*sech((z-80)/2.5).^2;

n2lab=n2lab*9.81;n2ocean=n2ocean*9.81;
n2ocean=(0.005/5)*sech((z-80)/5).^2+(0.0005/30)*sech((z-40)/30).^2;
subplot(1,2,1)
plot(rhobarlab,z,'k',rhobarocean,z,'k--')
grid on
legend('lab-like','ocean-like')
xlabel('density (kg/m^3)')
ylabel('z (m)')
subplot(1,2,2)
plot(n2lab,z,'k',n2ocean,z,'k--')
axis([-0.001 0.02 0 100])
grid on
xlabel('N^2 (1/s^2)')
