clear
close all

load IOxcyc.dat
load IOxcye.dat
load IOxeyc.dat
load surface1.dat

figure(1),
x1=IOxcyc(:,1);
y1=IOxcyc(:,2);
hold on
plot(surface1(:,1),surface1(:,2),'k-')
plot(x1,y1,'.r');
hold off
%axis([-0.6 0.6 -0.6 0.6])
title('xcyc');

figure(2),
x2=IOxcye(:,1);
y2=IOxcye(:,2);
hold on
plot(surface1(:,1),surface1(:,2),'k-')
plot(x2,y2,'.r');
hold off
%axis([-0.6 0.6 -0.6 0.6])
title('xcye');

figure(3),
x3=IOxeyc(:,1);
y3=IOxeyc(:,2);
hold on
plot(surface1(:,1),surface1(:,2),'k-')
plot(x3,y3,'.r');
hold off
%axis([-0.6 0.6 -0.6 0.6])
title('xeyc');