K = 500; % 100, 2000
mu0 = 4*pi*1e-7;
sigma = 6.99e6;
d = 1/3*1.5e-3;
omega = [10:10:200000]*2*pi;
%npts = 1000; omega = (1/npts:1/npts:1)*3125*2*pi;
%sigma = 2.17e6;
%d = 1.59e-3;
%omega = [25 60 200]*2*pi;

Delta = d*sqrt(K*mu0*omega*sigma/2);

Geddy = (1+1j)./Delta.*tan(Delta/(1+1j));
gain = abs(Geddy);
phase = angle(Geddy);


[g,p] = bode(tf(1,[1/2/pi/1800 1]),omega);

%figure;
ax(1)=subplot(211);
%semilogx(omega/2/pi, 20*log10([gain' squeeze(g)]));
semilogx(omega/2/pi, 20*log10([gain']));
hold all
ax(2)=subplot(212);
%semilogx(omega/2/pi,[phase'*180/pi squeeze(p)])
semilogx(omega/2/pi, [phase']*180/pi)
hold all
linkaxes(ax,'x')