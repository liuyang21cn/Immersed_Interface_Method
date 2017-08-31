clear
close all

load swo.dat
[ms,ns]=size(swo);
dtheta=360/(ns-1);
theta=0:dtheta:360;
figure,
hold on
for i=1:ms
plot(theta,swo(i,:));
end
hold off
xlabel('theta');ylabel('vorticity');
title('surface vorticity');
axis tight