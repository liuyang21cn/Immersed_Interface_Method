load gx.dat;
load gy.dat;
ip=2:5;

% figure(3)
% plot(gx(:,1),gx(:,ip),'-o')
% legend('u','v','p','o')
% 
% figure(4)
% plot(gy(:,1),gy(:,ip),'-o')
% legend('u','v','p','o')

load fbypp.dat;
cdcl=fbypp;
[m,n]=size(cdcl);
k=1;

figure(1)
plot(cdcl(1:k:m,1),cdcl(1:k:m,2),'r-')
xlabel('t')
ylabel('cd')

figure(2)
plot(cdcl(1:k:m,1),cdcl(1:k:m,3),'r-')
xlabel('t')
ylabel('cl')

