%% This code plots wave packets at a few times

% grid
L=130;
x=linspace(0,L,1001);
% wave paramaters
k=2*pi/3;
cp1=1.75
cp2=2.25
cg=2; cg1=cg,cg2=cg
t1=0
x0=20
t2=20
t3=40
packetlen=10
% functions
f11=exp(-((x-x0-cg1*t1)/packetlen).^2).*cos(k*(x-cp1*t1));
f12=exp(-((x-x0-cg1*t2)/packetlen).^2).*cos(k*(x-cp1*t2));
f13=exp(-((x-x0-cg1*t3)/packetlen).^2).*cos(k*(x-cp1*t3));
f21=exp(-((x-x0-cg2*t1)/packetlen).^2).*cos(k*(x-cp2*t1));
f22=exp(-((x-x0-cg2*t2)/packetlen).^2).*cos(k*(x-cp2*t2));
f23=exp(-((x-x0-cg2*t3)/packetlen).^2).*cos(k*(x-cp2*t3));



figure(1)
clf
betterplots
subplot(2,1,1)
plot(x,f11,'b',x,f12,'r',x,f13,'k')
hold on
plot([21 21],[-1 1],'bp--')
plot([21+cp1*(t2-t1) 21+cp1*(t2-t1)],[-1 1],'rp--')
plot([21+cp1*(t3-t1) 21+cp1*(t3-t1)],[-1 1],'kp--')
axis([0 130 -1 1])
grid on
ylabel('A')
text(3,0.8,'(a)')
subplot(2,1,2)
plot(x,f21,'b',x,f22,'r',x,f23,'k')
hold on
plot([21 21],[-1 1],'bp--')
plot([21+cp2*(t2-t1) 21+cp2*(t2-t1)],[-1 1],'rp--')
plot([21+cp2*(t3-t1) 21+cp2*(t3-t1)],[-1 1],'kp--')
axis([0 130 -1 1])
grid on
ylabel('A')
xlabel('x')
text(3,0.8,'(b)')


