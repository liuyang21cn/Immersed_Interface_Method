clear all

load xc.dat
load yc.dat
load p.dat

figure(1)
H=pcolor(xc,yc,p);
shading interp;
caxis([-1 1]);
axis equal;

figure(2)
surf(xc,yc,p)

