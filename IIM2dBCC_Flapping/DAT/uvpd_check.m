clear all
close all

load xc.dat
load yc.dat

load uc.dat
load vc.dat
load p.dat
load d.dat

load eu.dat
load ev.dat
load ep.dat

mc=size(xc,1);
nc=size(yc,1);

ixcs=1;
ixce=mc;
jycs=1;
jyce=nc;

figure(1)
surf(xc(ixcs:ixce),yc(jycs:jyce),uc(jycs:jyce,ixcs:ixce))
figure(2)
surf(xc(ixcs:ixce),yc(jycs:jyce),eu(jycs:jyce,ixcs:ixce))

figure(3)
surf(xc(ixcs:ixce),yc(jycs:jyce),vc(jycs:jyce,ixcs:ixce))
figure(4)
surf(xc(ixcs:ixce),yc(jycs:jyce),ev(jycs:jyce,ixcs:ixce))

figure(5)
surf(xc(ixcs:ixce),yc(jycs:jyce),p(jycs:jyce,ixcs:ixce))
figure(6)
surf(xc(ixcs:ixce),yc(jycs:jyce),ep(jycs:jyce,ixcs:ixce))

figure(7)
surf(xc(ixcs:ixce),yc(jycs:jyce),d(jycs:jyce,ixcs:ixce))

