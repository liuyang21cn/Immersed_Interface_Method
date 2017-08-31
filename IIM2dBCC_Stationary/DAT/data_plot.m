clear all

load xc.dat
load yc.dat

ms=2;
for k=1:ms
  if k==1
     load surface1.dat
  elseif k==2
     load surface2.dat
  elseif k==3
     load surface3.dat
  end
end

kip=1;

icontour=1;
iplot=1;
ishape=1;
iquiver=0;
iplotall=0;

if iplotall==1
   load uc.dat
   load vc.dat
   load p.dat
   load d.dat
end

if iquiver==1
  load uc.dat
  load vc.dat
end

if iplot==1
  load wo.dat
elseif iplot==2
  load ph.dat
elseif iplot==3
  load p.dat
elseif iplot==4
  load d.dat
end

igiven=0;
label=0;
level=100;

mc=size(xc,1);
nc=size(yc,1);

ixcs=1;
ixce=mc;
jycs=1;
jyce=nc;

wov=[-2.5:0.1:2.5];
phv1=[-8:0.2:8];
phv2=[-0.02:0.002:0.02];
phv=[phv1,phv2];
pv=[-4:0.2:4];

figure(1)
if icontour==1
  if iplot==1
    if igiven==0
      H=pcolor(xc,yc,wo);
      shading interp;
      caxis([-5 5]);
      axis equal;
    else
      cs=contour(xc,yc,wo,wov);
    end
  end
  if iplot==2
    if igiven==0
      cs=contour(xc,yc,ph,level);
    else
      cs=contour(xc,yc,ph,phv);
    end
  end
  if iplot==3
    if igiven==0
      H=pcolor(xc,yc,p);
      shading interp;
      caxis([-1 1]);
      axis equal;
%      [cs,h]=contour(xc,yc,p,level);
%      surf(xc,yc,p)
%      set(h,'ShowText','on','TextStep',get(h,'LevelStep')*4)
    else
      cs=contour(xc,yc,p,pv);
    end
  end  
  if iplot==4
%    H=pcolor(xc,yc,d);
%    shading interp;
%    caxis([-1 1]);
%    axis equal;
    cs=contour(xc,yc,d,level);
  end
  if label==1
    clabel(cs)
  end
  hold on
end  
if iquiver==1
  quiver(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),...
         uc(jycs:kip:jyce,ixcs:kip:ixce),vc(jycs:kip:jyce,ixcs:kip:ixce))
  hold on
end
if ishape==1
  for k=1:ms
    if k==1
      plot(surface1(:,1),surface1(:,2),'k-')
      hold on
    elseif k==2
      plot(surface2(:,1),surface2(:,2),'k-')
      hold on
    elseif k==3
      plot(surface3(:,1),surface3(:,2),'k-')
    end
  end
end
hold off
axis equal
% axis([-1 2 -1 1])
xlabel('x')
ylabel('y')

if iplotall==1

  figure(2)
  surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),uc(jycs:kip:jyce,ixcs:kip:ixce))

  figure(3)
  surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),vc(jycs:kip:jyce,ixcs:kip:ixce))

  figure(4)
  surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),p(jycs:kip:jyce,ixcs:kip:ixce))

  figure(5)
  surf(xc(ixcs:kip:ixce),yc(jycs:kip:jyce),d(jycs:kip:jyce,ixcs:kip:ixce))

end

