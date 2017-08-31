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
hold on
plot(cdcl(1:k:m,1),cdcl(1:k:m,2),'r-');
plot(cdcl(1:k:m,1),cdcl(1:k:m,4),'b-')
xlabel('t')
ylabel('cd')

figure(2)
hold on
plot(cdcl(1:k:m,1),cdcl(1:k:m,3),'r-')
plot(cdcl(1:k:m,1),cdcl(1:k:m,5),'b-')
xlabel('t')
ylabel('cl')

