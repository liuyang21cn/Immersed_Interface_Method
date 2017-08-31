clear all

load xc.dat
load yc.dat

load gx.dat
load gy.dat

load p.dat
load surface1.dat

kip=1;
level=60;

mc=size(xc,1);
nc=size(yc,1);

ixcs=1;
ixce=mc;
jycs=1;
jyce=nc;

figure(1)
plot(gx(:,1),gx(:,4),'rx-',gy(:,1),gy(:,4),'bo-')

figure(2)
contour(xc,yc,p,level)
hold on
plot(surface1(:,1),surface1(:,2),'k-')
axis equal
hold off


